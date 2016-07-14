#!/usr/bin/env python

from __future__ import print_function

__author__ = 'konradjk'

import argparse
import subprocess
import sys
import gzip
import os
import os.path
import time
import re
import pipes
import tempfile

# Great hack for 2.X and 3.X to use input()
try:
    input = raw_input
except NameError:
    pass


def main(args, pass_through_args):
    if args.lookup:
        if args.dbnsfp or args.gtex or args.ancestral or args.context or args.add_rank_filter \
                or args.filter is not None or args.lof_only or args.refseq or args.basic or args.apply_all \
                or args.skip_conservation or args.position != 0.05 or args.intron_size != 15:
            print('Note that --lookup can only be used with the default options. Please re-run without this option.', file=sys.stderr)
            sys.exit(1)
    bhosts = subprocess.check_output(["bhosts"], stderr=subprocess.STDOUT)

    if args.cache_version != args.vep_version:
        print('WARNING: Cache version is not the same as VEP version. Continuing, but results may be off...', file=sys.stderr)
    if '.vcf' not in args.vcf:
        print('ERROR: ".vcf" is not in your input file - is this a VCF?', file=sys.stderr)
        sys.exit(1)
    fasta_type = ''
    try:
        subprocess.check_output(['samtools'], stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        if 'Version: 0.1' in e.output:
            print('Samtools version 0.1 detected.', file=sys.stderr)
            fasta_type = '.rz'
        elif 'Version: 1' in e.output:
            print('Samtools version 1 detected.', file=sys.stderr)
            fasta_type = '.gz'
    except OSError:
        print('ERROR: You do not have samtools on your path. Needed for LOFTEE! (Run: use Samtools, and/or put it in your ~/.my.bashrc)', file=sys.stderr)
        sys.exit(1)
    try:
        subprocess.check_output(['git'], stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        pass
    except OSError:
        print('ERROR: You do not have git on your path. Please run: use Git-2.0, and/or put it in your ~/.my.bashrc')
        sys.exit(1)
    if args.mysql or args.lookup:
        conf_file = os.path.join(os.path.expanduser('~'), '.my.cnf')
        if not os.path.exists(conf_file) or '[loftee]' not in open(conf_file).read():
            print('ERROR: [loftee] tag not found in ~/.my.cnf - please add this tag', file=sys.stderr)
            sys.exit(1)

    konrad_vep_dir_raw = '/humgen/atgu1/fs03/konradk/vep/%s'
    loftee_data_location = args.loftee_data_location
    combine_script_path = os.path.join(os.path.dirname(__file__), "combine_and_delete.py")
    check_script_path = os.path.join(os.path.dirname(__file__), "check_vep_jobs.py")
    vep_location = args.vep_location or '/humgen/atgu1/fs03/konradk/vep/ensembl-tools-release-%s/scripts/variant_effect_predictor/variant_effect_predictor.pl' % args.vep_version
    vep_cache_dir = args.vep_cache_dir if args.vep_cache_dir else konrad_vep_dir_raw % ''
    dev = '_dev' if args.dev else ''
    git_version = subprocess.check_output(['git', '--git-dir=%s/loftee%s/.git/' % (konrad_vep_dir_raw % '', dev), 'describe', '--tags']).strip().split('-')[0]
    num_known = 0
    if args.no_split:
        file_no = 1
        out_base = os.path.basename(args.vcf)
        all_current_line = 'all'
    else:
        cat_file = 'zcat' if args.vcf.endswith('.gz') else 'cat'

        lines = None
        # Make decision on how many lines per file
        if args.split_size is not None:
            number_to_split = args.split_size if args.split_size > 0 else 3E9
        elif args.default_submission:
            number_to_split = 20000
        else:
            print("Getting number of lines...", file=sys.stderr)
            p1 = subprocess.Popen([cat_file, args.vcf], stdout=subprocess.PIPE)
            p2 = subprocess.Popen(["wc", "-l"], stdin=p1.stdout, stdout=subprocess.PIPE)
            p1.stdout.close()
            lines = int(p2.communicate()[0])

            print("File is %s lines." % lines, file=sys.stderr)
            number_to_split = 20000
            print("Typical safe limits: <=10K variants on the hour queue; <=500K variants on the priority queue", file=sys.stderr)
            try:
                s = input('File is %s lines. How many lines per split file? [%s] ' % (lines, number_to_split))
                number_to_split = int(s)
                if number_to_split < 1:
                    number_to_split = lines + 1
            except Exception as e:
                pass

        # Prepare files
        out_base = os.path.basename(args.vcf)
        file_no = 0
        current_line = 0
        all_current_line = 0
        header = ''
        my_open = gzip.open if args.vcf.endswith('.gz') else open
        g = None
        known = None
        header_line = None
        try:
            os.makedirs(args.output)
        except Exception as e:
            pass

        if lines is not None:
            print("Splitting into ~%s files of %s lines..." % (lines/number_to_split + 1, number_to_split), file=sys.stderr)
        else:
            print("Splitting into files of %s lines..." % (number_to_split), file=sys.stderr)

        start = time.time()
        if args.lookup:
            import MySQLdb
            import MySQLdb.cursors

            db = MySQLdb.connect(read_default_group='loftee', cursorclass=MySQLdb.cursors.DictCursor)
            conn = db.cursor()
        with my_open(args.vcf) as f:
            for line in f:
                if line.startswith('#'):
                    if line.startswith('#CHROM'):
                        header_line = line.lstrip('#').strip().split('\t')
                        header_line = dict(zip(header_line, range(len(header_line))))
                    header += line
                else:
                    if args.lookup:
                        if known is None:
                            if header_line is None:
                                print("Header line not found in VCF and need it to lookup entries", file=sys.stderr)
                                sys.exit(1)

                            pipe = pipes.Template()
                            pipe.append('bgzip -c /dev/stdin', '--')
                            base = out_base.replace('.vcf', '_known.vep.vcf')
                            if not base.endswith('.gz'): base += '.gz'
                            known = pipe.open(os.path.join(args.output, base), 'w')
                            known.write(header)
                        fields = line.strip().split('\t')
                        sql_data = [args.vep_version, git_version[:10], fields[header_line['CHROM']], fields[header_line['POS']], fields[header_line['REF']].upper(), fields[header_line['ALT']].upper()]
                        conn.execute('SELECT csq from known_variants WHERE vep=%s AND loftee=%s AND chrom=%s AND pos=%s AND ref=%s AND alt=%s', sql_data)
                        data = conn.fetchone()
                        if data is not None:
                            info = ';'.join([x for x in re.split(';(?=\w)', fields[header_line['INFO']].rstrip(';')) if x.split('=')[0] != 'CSQ'])
                            fields[header_line['INFO']] = info + ';CSQ=' + data['csq']
                            known.write('\t'.join(fields) + '\n')
                            num_known += 1
                            continue
                    if g is None:
                        base = out_base.replace('.vcf', '_%04d.vcf' % file_no)
                        if not base.endswith('.gz'): base += '.gz'
                        g = gzip.open(os.path.join(args.output, base), 'w')
                        file_no += 1
                        g.write(header)
                    if args.sites_only:
                        g.write(line.strip().split('\t')[:8])
                    else:
                        g.write(line)
                    current_line += 1
                    all_current_line += 1
                    if current_line == number_to_split:
                        g.close()
                        g = None
                        current_line = 0
        if g is not None: g.close()
        if known is not None: known.close()
        if not file_no:
            base = out_base.replace('.vcf', '_known.vep.vcf')
            if not base.endswith('.gz'): base += '.gz'
            known_file = os.path.join(args.output, base)
            final_name = known_file.replace('_known', '')
            os.rename(known_file, final_name)
            subprocess.check_output(['tabix', final_name])
            print('Yay! All your variants have been seen before. Annotated file is now ready at: %s' % os.path.abspath(final_name))
            sys.exit(0)

        print("Done. Wrote %s files of %s lines (took %s seconds)." % (file_no, number_to_split, time.time() - start), file=sys.stderr)

    possible_queues = {'w': 'week',
                       'p': 'priority',
                       'h': 'hour',
                       'g': 'gsa'}
    queue = 'hour'
    if args.queue is None:
        if args.default_submission:
            args.queue = 'hour'
        else:
            args.queue = input('Which queue to submit to: (w)eek, (p)riority, (h)our? [hour] ')

    if len(args.queue.strip()) > 0 and args.queue[0] in possible_queues:
        queue = possible_queues[args.queue[0]]

    additional_options = ',filter_position:%s' % float(args.position)
    additional_options += ',min_intron_size:%s' % int(args.intron_size)
    if args.apply_all: additional_options += ',apply_all:true'

    if not args.skip_conservation:
        if args.mysql:
            additional_options += ',conservation_file:mysql'
        else:
            additional_options += ',conservation_file:%s/phylocsf.sql' % loftee_data_location

    log_dir = os.path.abspath(os.path.join(args.output, 'logs'))
    try:
        os.makedirs(log_dir)
    except Exception as e:
        pass

    # Set environment vars
    env = os.environ.copy()
    loftee_dir = konrad_vep_dir_raw % ('loftee_dev/' if args.dev else 'loftee/',)
    env['PERL5LIB'] = ':'.join([loftee_dir, env['PERL5LIB']]) if 'PERL5LIB' in env else loftee_dir
    print("PERL5LIB will be: %s" % env['PERL5LIB'], file=sys.stderr)

    print("Submitting %s jobs to the %s queue..." % (file_no, queue), file=sys.stderr)
    started_time = time.strftime("%Y_%m_%d_%H_%M_%S")
    log = open(os.path.join(log_dir, 'starting_log_%s.txt' % started_time), 'w')
    log.write('#command_line_call=%s\n' % ' '.join(sys.argv))
    log.write('#current_working_directory=%s\n' % os.getcwd())
    log.write('#loftee_version=%s\n' % git_version)
    log.write('#vep_version=%s\n' % args.vep_version)
    log.write('#cache_version=%s\n' % args.cache_version)
    log.write('#original_file=%s\n' % os.path.abspath(args.vcf))
    log.write('#time_started=%s\n' % started_time)
    log.write('#original_size=%s\n' % all_current_line)
    log.write('#number_of_files=%s\n' % file_no)
    if args.lookup: log.write('#number_of_known_variants=%s\n' % num_known)

    project = os.path.basename(args.vcf) if args.project is None else args.project

    start = time.time()
    all_jobs = []
    error_log = None
    for i in range(file_no):
        if args.no_split:
            input_file = os.path.abspath(args.vcf)
        else:
            base = out_base.replace('.vcf', '_%04d.vcf' % i)
            if not base.endswith('.gz'): base += '.gz'
            input_file = os.path.abspath(os.path.join(args.output, base))
        if file_no > 1 or args.lookup:
            base = out_base.replace('.vcf', '_%04d.vep.vcf' % i)
        else:
            base = out_base.replace('.vcf', '.vep.vcf')
        if not base.endswith('.gz'): base += '.gz'
        output_file = os.path.abspath(os.path.join(args.output, base))

        cache_dir_sub = args.cache_version if args.cache_version <= 75 else '%s_GRCh%s' % (args.cache_version, args.assembly_version)
        fasta_sub = args.assembly_version if args.assembly_version > 37 else '%s.%s' % (args.assembly_version, min(args.cache_version, 75))
        fasta_location = '%s/homo_sapiens/%s/Homo_sapiens.GRCh%s.dna.primary_assembly.fa' % (vep_cache_dir, cache_dir_sub, fasta_sub)

        job_name = out_base.replace('.vcf', '_%04d' % i).replace('.gz', '')
        bsub_command = ['bsub', '-oo', os.path.join(log_dir, 'log_%04d.txt' % i), '-R', 'rusage[mem=%s]' % args.memory, '-q', queue]
        temp_directory = tempfile.mkdtemp()
        if args.vep_version >= 79:
            bsub_command.extend(['-E', 'ls %(cache_dir)s > /dev/null; mkdir -p %(temp)s; ln -f -s %(fasta)s %(temp)s' % {'cache_dir': vep_cache_dir, 'fasta': fasta_location, 'temp': temp_directory}])
        else:
            bsub_command.extend(['-E', 'cd %s' % os.getcwd()])
        bsub_command.extend(['-g', '/%s' % project, '-P', project, '-J', job_name])
        if file_no == 1 and not args.no_email and not args.lookup:
            bsub_command.append('-N')
        if queue == 'hour':
            bsub_command.extend(['-W', '4:00'])

        if args.high_io is not None: bsub_command.extend(['-R', 'rusage[%s]' % args.high_io])

        bsub_command.append('/broad/software/free/Linux/redhat_5_x86_64/pkgs/perl_5.10.1/bin/perl')
        bsub_command.append(vep_location)

        # VEP options
        if not args.skip_everything: bsub_command.append('--everything')
        bsub_command.extend(['--vcf', '--allele_number', '--no_stats'])
        bsub_command.extend(['--cache', '--offline', '--dir', vep_cache_dir, '--force_overwrite'])
        bsub_command.extend(['--cache_version', str(args.cache_version)])
        if args.force_vcf: bsub_command.extend(['--format', 'vcf'])
        if args.filter: bsub_command.extend(['--filter', args.filter])

        if args.vep_version >= 79:
            bsub_command.extend(['--fasta', '%s/%s' % (temp_directory, os.path.basename(fasta_location))])
        else:
            bsub_command.extend(['--fasta', fasta_location])
        if args.vep_version >= 80: bsub_command.append('--minimal')
        if args.cache_version > 75: bsub_command.extend(['--assembly', 'GRCh%s' % args.assembly_version])

        if output_file.endswith('.gz'):
            bsub_command.append('--tabix')

        if args.basic: bsub_command.append('--gencode_basic')
        if args.refseq: bsub_command.append('--refseq')

        # Plugins
        if not args.skip_lof: bsub_command.extend(['--plugin', 'LoF,human_ancestor_fa:%s/human_ancestor.fa%s%s' % (loftee_data_location, fasta_type, additional_options)])

        if args.lof_only:
            bsub_command.extend(['--plugin', 'RankFilter,initiator_codon_variant'])
        elif args.add_rank_filter:
            bsub_command.extend(['--plugin', 'RankFilter,intron_variant'])
        if args.gtex: bsub_command.extend(['--plugin', 'TissueExpression,db_location:%s/gtex.db' % loftee_data_location])
        if args.dbnsfp is not None: bsub_command.extend(['--plugin', 'dbNSFP,%s,%s' % (args.dbnsfp_data_location, args.dbnsfp)])
        if args.context is not None: bsub_command.extend(['--plugin', 'context,%s' % args.context])
        if args.ancestral: bsub_command.extend(['--plugin', 'ancestral,human_ancestor_fa:%s/human_ancestor.fa%s' % (loftee_data_location, fasta_type)])
        if args.dbscSNV: bsub_command.extend(['--plugin', 'dbscSNV,%s/misc/dbscSNV.txt.gz' % loftee_data_location])
        bsub_command.extend(pass_through_args)
        bsub_command.extend(['-i', input_file])
        bsub_command.extend(['-o', output_file])


        print("Running: " + " ".join(bsub_command), file=sys.stderr)

        if not args.dry:
            job_number = run_job_lsf(bsub_command, env)
            if job_number is None:
                if error_log is None:
                    error_log = open(os.path.join(log_dir, 'error_log_%s.txt' % started_time), 'w')
                print(' '.join(['"%s"' % x if ' ' in x else x for x in bsub_command]), file=error_log)
            else:
                log.write('%s\t%s' % (job_number, "\t".join(map(str, bsub_command)) + '\n'))
                all_jobs.append(job_number)

    if file_no > 1 or args.lookup:
        base = os.path.abspath(os.path.join(args.output, out_base.replace('.vcf', '_%04d.vep.vcf')))
        output = os.path.abspath(os.path.join(args.output, out_base.replace('.vcf', '.vep.vcf')))
        if not base.endswith('.gz'):
            base += '.gz'
            output += '.gz'

        bsub_command = ['bsub', '-oo', os.path.join(log_dir, 'log_done.txt')]
        if args.lookup:
            bsub_command.extend(['-q', queue if queue == 'priority' else 'week'])
        else:
            bsub_command.extend(['-q', queue])
            if queue == 'hour': bsub_command.extend(['-W', '4:00'])
        bsub_command.extend(['-E', 'python %s --base %s --number %s' % (check_script_path, base, file_no), '-g', '/%s' % project, '-P', project])
        if not args.no_email: bsub_command.append('-N')
        bsub_command.extend(['-w', ' && '.join(all_jobs[-1000:])]) # limiting to last 1K jobs since bash doesn't like long argument lists
        bsub_command.extend(['python', combine_script_path])

        bsub_command.extend(['--base', base, '--output', output])
        if not args.no_delete: bsub_command.append('--delete')
        if args.lookup: bsub_command.extend(['--lookup', '%s,%s' % (args.vep_version, git_version[:10])])
        if args.skip_lookup_add: bsub_command.append('--skip_lookup_add')
        bsub_command.extend(['--number', str(file_no)])

        print("Running: " + " ".join(bsub_command), file=sys.stderr)

        if not args.dry:
            job_number = run_job_lsf(bsub_command, env)
            log.write('%s\t%s\n' % (job_number, "\t".join(map(str, bsub_command))))

            # Check if any jobs have failed
            new_queue = 'hour' if queue != 'priority' else queue
            dependency = 'ended(%s) && (exit(%s,>0))' % (') && ended('.join(all_jobs), ',>0) || exit('.join(all_jobs))
            bsub_command = ['bsub', '-N', '-q', new_queue, '-w', dependency, '-g', '/%s' % project, '-P', project, '-J', '%s_chk' % project, 'echo']
            bsub_command.append('"At least one job failed. Please rerun with: cd %s; python /humgen/atgu1/fs03/konradk/src/rerun_failed_jobs.py -o %s [--dry]"' % (os.getcwd(), args.output))
            rerun_job_number = run_job_lsf(bsub_command, env)

            # If final job succeeded, kill fail checkpoint
            bsub_command = ['bsub', '-o', '/dev/null', '-q', new_queue, '-g', '/%s' % project, '-P', project, '-J', '%s_chk2' % project, '-w', job_number, 'bkill', rerun_job_number]
            run_job_lsf(bsub_command, env)

    if args.lookup:
        conn.close()
        db.close()
    log.close()
    if error_log is not None: error_log.close()

    print("Done submitting %s jobs! Took %s seconds." % (len(list(range(file_no))) + 1, time.time() - start), file=sys.stderr)


def run_job_lsf(bsub_command, env=None):
    try:
        bsub_output = subprocess.check_output(bsub_command, stderr=subprocess.STDOUT, env=env)
        return [x for x in bsub_output.split('\n') if x.startswith('Job')][0].split('<')[1].split('>')[0]
    except subprocess.CalledProcessError as e:
        print(e, file=sys.stderr)
        print(e.output, file=sys.stderr)
        print(' '.join(['"%s"' % x if ' ' in x else x for x in bsub_command]), file=sys.stderr)
        return None

if __name__ == '__main__':
    INFO = """Runs VEP with LOFTEE.
Minimal usage is: python run_lof_annotation.py -i input.vcf[.gz] -o output_directory
This will count the number of lines in the file, and prompt for the number of lines to split the file into and which queue to submit to.
If you'd like to omit the interactive step, add -s NUMBER_OF_LINES -q QUEUE. -s 20000 -q hour is a sensible default."""

    try:
        import configargparse
        parser = configargparse.ArgumentParser(description=INFO,
            default_config_files=['~/.run_lof_config'],
            args_for_setting_config_path=["-c", "--run-lof-config"],
            formatter_class=configargparse.DefaultsRawFormatter)
    except ImportError:
        parser = argparse.ArgumentParser(description=INFO, formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument("-c", "--run-lof-config", dest="config", help="To enable this option, please install configargparse")

    parser.add_argument('--vcf', '--input', '-i', help='Input VCF file; may be gzipped', required=True)
    parser.add_argument('--output', '-o', help='Output directory', required=True)

    submission_options = parser.add_argument_group('Submission options')
    submission_options.add_argument('--split_size', '-s', help='Number of lines to split file into', type=int)
    submission_options.add_argument('--queue', '-q', help='Which queue to submit to: (h)our, (p)riority, (w)eek')
    submission_options.add_argument('--default_submission', '-d', help='Submit default options (20000 variants to hour queue)', action='store_true')
    submission_options.add_argument('--memory', help='GB of memory to submit', type=int, default=8)
    submission_options.add_argument('--high_io', help='Set I/O requirements (default none, if used as flag, will default to /humgen/atgu1/fs03/ high I/O)', nargs='?', const='argon_io=1000', default=None)
    submission_options.add_argument('--no_email', action='store_true', help='Do not send an email when done')
    submission_options.add_argument('--no_delete', action='store_true', help='Do not delete intermediate files')
    submission_options.add_argument('--no_split', "-n", help='Do not split file, but run VEP/LOFTEE on one file', action='store_true')
    submission_options.add_argument('--project', '-P', help='Project name for submission [default = input vcf name]')
    submission_options.add_argument('--dry', help='Dry run (creates directories and splits file, but does not submit any jobs)', action='store_true')
    submission_options.add_argument('--server', help="Deprecated. This option is ignored.")
    submission_options.add_argument('--lookup', help="Use annotation lookup to speed up", action='store_true')
    submission_options.add_argument('--skip_lookup_add', help="Do not add new annotations to the database.", action='store_true')
    submission_options.add_argument('--sites_only', help="Create sites VCFs for annotation.", action='store_true')

    vep_options = parser.add_argument_group('VEP options')
    vep_options.add_argument('--cache_version', help='Default: 79', default=79, type=int)
    vep_options.add_argument('--vep_version', help='Default: 79', default=79, type=int)
    vep_options.add_argument('--assembly_version', help='Default: 37 (i.e. GRCh37); can also be 38', default=37, type=int)
    vep_options.add_argument('--force_vcf', help='Force VEP to recognize input file as a VCF', action='store_true')
    vep_options.add_argument('--basic', '--gencode_basic', help='Only use transcripts from Gencode Basic annotation', action='store_true')
    vep_options.add_argument('--refseq', help='Use RefSeq transcript models for annotation', action='store_true')
    vep_options.add_argument('--lof_only', help='Only print possible LoF (initiator_codon_variant and above) variants in output VCF', action='store_true')
    vep_options.add_argument('--skip_everything', help='Skip --everything flag to VEP', action='store_true')
    vep_options.add_argument('--filter', help='Comma-separated list of consequences (eg. SO terms) to pass to the VEP '
        '--filter option. The CSQ key/value will only be added to the INFO field for varaints that *do* have one of '
        'these consequences. For a list of valid consequence terms, see '
        'http://useast.ensembl.org/info/genome/variation/predicted_data.html#consequences')
    additional_plugins = parser.add_argument_group('Additional plugins')
    additional_plugins.add_argument('--add_rank_filter', help='Filter for intron variant and up (default: run VEP with --plugin RankFilter,intron_variant)', action='store_true')
    additional_plugins.add_argument('--context', help='Run context plugin', type=int, nargs='?', const=1, default=None)
    additional_plugins.add_argument('--ancestral', help='Run ancestral plugin', action='store_true')
    additional_plugins.add_argument('--gtex', help='Add GTEx information', action='store_true')
    additional_plugins.add_argument('--dbnsfp', help='Add dbNSFP annotations. The value should be a comma-separated list of dbNSFP column names or nothing (eg. use as a flag) to only print all available dbNSFP columns.', nargs='?', const="")
    additional_plugins.add_argument('--dbscSNV', '--dbscsnv', help='Add dbscSNV annotations', action='store_true')

    loftee_options = parser.add_argument_group('LOFTEE options')
    loftee_options.add_argument('--skip_lof', help='Skip LOFTEE', action='store_true')
    loftee_options.add_argument('--position', help='Position in transcript where a variant should be filtered (filter_position in LOFTEE)', type=float, default=0.05)
    loftee_options.add_argument('--intron_size', help='Minimum intron size, below which a variant should be filtered (min_intron_size in LOFTEE)', type=int, default=15)
    loftee_options.add_argument('--skip_conservation', help='Skip PhyloCSF filter', action='store_true')
    loftee_options.add_argument('--apply_all', help='Apply LoF filters to all variants, not just LoF', action='store_true')
    loftee_options.add_argument('--mysql', help='Uses MySQL instead of SQLite for PhyloCSF', action='store_true')
    loftee_options.add_argument('--dev', help='Runs development version of LOFTEE (not recommended)', action='store_true')

    location_options = parser.add_argument_group('File paths')
    location_options.add_argument("--vep_location", help="Full path of variant_effect_predictor.pl")
    location_options.add_argument("--loftee_data_location", default='/humgen/atgu1/fs03/konradk/loftee_data/')
    location_options.add_argument("--vep_cache_dir", help="VEP cache directory")
    location_options.add_argument("--dbnsfp_data_location", default='/humgen/atgu1/fs03/weisburd/xbrowse/data/reference_data/dbNSFP/dbNSFPv2.9.gz')

    args, pass_through_args = parser.parse_known_args()

    if parser.__module__ == "configargparse":
        print("Running with the following settings and default values: ")
        parser.print_values()
    elif args.config:
        parser.error("To enable the -c / --run-lof-config option, please run: pip install --user configargparse")
    main(args, pass_through_args)

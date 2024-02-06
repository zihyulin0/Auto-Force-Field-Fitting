# Author: Brett Savoie and Zih-Yu Lin
import os, sys, subprocess, argparse
import writing_template as wt
import common_static as cs

FLAG_COMMAND = 'command'
FLAG_JOBNAME = 'jobname'
FLAG_PROCS = '-procs'
FLAG_WT = '-walltime'
FLAG_QUEUE = '-queue'
FLAG_PPN = '-ppn'

class ShellSubmit:
    """
    submit a shell command  (like a python call) to slurm
    """
    def __init__(self, options):
        """
        shell submit

        :type options: parser object
        :param options: all user formatted input
        """
        self.options = options
        self.options.outname = '{}.submit'.format(self.options.jobname)
        self.working_dir = os.getcwd()
        self.jobid = None

    def WriteSubmission(self):
        """
        Write slurm sbatch script
        """
        # due to duplicate method names
        header =  wt.Writing_Template(self.options)

        # Write the header
        header.WriteSbatchHeader()
        # Write custom commands
        with open(self.options.outname,'a') as f:
            f.write("\n# cd into the submission directory\n")
            f.write("cd -- {}\n".format(self.working_dir))
            f.write("echo Working directory is {}\n".format(self.working_dir))

            f.write("# Copy the local path\n\n")

            f.write("# Run python script\n")
            f.write("sleep 2\n")
            f.write("{}\n".format(self.options.command))
            f.write("sleep 2\n")

    def run(self):
        """
        run the workflow
        """
        self.WriteSubmission()
        self.jobid = cs.SubmitGetJobID(self.options.outname)

def main(args):
    # get, and format inputs
    options = parse_and_format_input(args)

    # Run the workflow
    driver =  ShellSubmit(options)
    driver.run()
    print('Job {} submitted'.format(driver.jobid))

    return driver.jobid


def get_parser():
    """
    Argparse wrappers
    Define user input options here

    :rtype: argparser object
    :return parser
    """

    parser = argparse.ArgumentParser(
        description='Submit a shell command to server')

    #positional arguments
    parser.add_argument(
        FLAG_COMMAND,
        default = '',
        help=' command to execute such a python foo.py -b 1',
        type=str)

    parser.add_argument(
        FLAG_JOBNAME,
        help=' prefix for sh and error name',
        type=str)

    # optional arguments
    parser.add_argument(
        FLAG_PROCS,
        default=1,
        help='proces to run',
        type=int)

    # this needs to be strings to toggle minute flag
    parser.add_argument(
        FLAG_WT,
        default='4',
        help='wall time',
        type=str)

    parser.add_argument(
        FLAG_QUEUE,
        default='bsavoie',
        help='queue',
        type=str)

    parser.add_argument(
        FLAG_PPN,
        default=1,
        help='ppn',
        type=int)

    return parser

def parse_and_format_input(args):
    """
    Format input into useful type

    :type options: argparse parser object
    :param options: object with all the inputs

    :rtype: parser object
    :return: formmated parser object
    """

    parser = get_parser()
    options = parser.parse_args(args)

    # deal with minute
    if 'min' in options.walltime:
        options.min_flag = 1
        options.walltime = int(options.walltime.strip('min'))
    else:
        options.min_flag = 0
        options.walltime = int(options.walltime)

    return options



if __name__ == '__main__':
    main(sys.argv[1:])

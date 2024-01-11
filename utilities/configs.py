import os, sys, re, getpass, shutil
from types import SimpleNamespace
from Excpetions import TAFFIException
from utilities.parse import read_n_column

class IOConfigException(TAFFIException):
    pass
class InputConfig:
    # configuration class that holds various params
    def __init__(self):
        self._running_dir = os.getcwd()
        # environment options
        self.env = self.InitEnvDefault()
        # force field fitting options
        self.ff_fit = self.InitFitConfigDefault()
        # geometry optimizatioin
        self.param_geoopt = SubmitJobConfig()
        # for bonds and angles
        self.param_ba = SubmitJobConfig()
        # for dihedrals
        self.param_d = SubmitJobConfig()
        # fitting
        self.fit = SubmitJobConfig()
        # charges md
        self.charges_md = SubmitJobConfig()
        # charges dft
        self.charges_qc = SubmitJobConfig()
        self.vdw_md = SubmitJobConfig()
        self.vdw_qc = SubmitJobConfig()

    def InitEnvDefault(self):
        # initialize enviroment default options
        return SimpleNamespace(
            TAFFI_PATH='/'.join(os.path.abspath(__file__).split('/')[:-2]), # TAFFI folder path
            USER = getpass.getuser(),
            RUN_DIR = self._running_dir,
            # different version of lammps might be problematic with mpi
            LAMMPS_EXE = '/depot/bsavoie/apps/lammps/exe/lmp_mpi_180501', # lammps executable
            # Use 4.1.2
            ORCA_EXE = '/depot/bsavoie/apps/orca_4_1_2/orca', # orca executable
            MODULE_STRING = '')

    def InitFitConfigDefault(self):
        # initialize force field fitting options
        return SimpleNamespace(
            # Create an empty db in the running folder if no existing db is
            # provided
            FF= '{}/empty.db'.format(self._running_dir), # existing FF db file
            CHARGE= 0,
            GENS= 2,
            BASIS= 'def2-TZVP',
            FUNCTIONAL= 'wB97X-D3',
            D3STR= '',
            XYZ = [f for f in os.listdir('.') if os.path.isfile(f) and "xyz" in f])

    def Parse(self, infile):
        PARAM_DB = '{}/Params_for_Batch.db'.format(self._running_dir)

        if infile != '' and not os.path.isfile(infile):
            raise IOConfigException('Cannot find specified config in file: {}'.format(infile))

        # return default configuration if no input configuration
        if infile == '':
            print('No configuration file specified, using default parameters...')
            with open(self.ff_fit.FF, 'w') as f:
                f.write('#emtpy db\n')
            shutil.copyfile(self.ff_fit.FF, PARAM_DB)
            self.ff_fit.FF = "{}/{}".format(self._running_dir, PARAM_DB)
            return

        data = read_n_column(infile, n=2)
        # key: option name
        # all in str
        dict_from_data = { row[0]: row[1] for row in data}

        env_attr = list(self.env.__dict__.keys())
        ff_fit_attr = list(self.ff_fit.__dict__.keys())
        job_attr = list(self.charges_md.__dict__.keys())
        prefixes = ['param_geoopt', 'param_ba', 'param_d', 'fit', 'charges_md', 'charges_qc', 'vdw_md', 'vdw_qc']

        for key, value in dict_from_data.items():
            if key in env_attr:
                setattr(self.env, key, value)
            elif key in ff_fit_attr:
                setattr(self.env, key, value)
            elif key.split('_')[-1] in job_attr:
                fields = key.split('_')
                prefix = '_'.join(fields[:-1]).lower()
                sub_attr = fields[-1]
                if prefix in prefixes:
                    setattr(getattr(self, prefix), sub_attr, value)

    def WriteConfig(self, outfile):
        format_str = "{:<20s} {:}\n"
        with open(outfile,'w') as f:
            for key, value in self.env.__dict__.items():
                f.write(format_str.format(key,value))
            for key, value in self.ff_fit.__dict__.items():
                f.write(format_str.format(key,value))
            prefixes = ['param_geoopt', 'param_ba', 'param_d', 'fit', 'charges_md', 'charges_qc', 'vdw_md', 'vdw_qc']
            for prefix in prefixes:
                for key, value in getattr(self, prefix).__dict__.items():
                    f.write(format_str.format(prefix.upper()+'_'+key, value))
    def Validate(self):
        """
        validate the options read in and check if the exes exist and deal with Param.db checks
        """
        PARAM_DB = '{}/Params_for_Batch.db'.format(self._running_dir)

        # fixing types
        prefixes = ['param_geoopt', 'param_ba', 'param_d', 'fit', 'charges_md', 'charges_qc', 'vdw_md', 'vdw_qc']
        # change str to int for job attributes
        for prefix in prefixes:
            getattr(self, prefix).FixType()
        self.ff_fit.CHARGE = float(self.ff_fit.CHARGE)
        self.ff_fit.GENS = int(self.ff_fit.GENS)

        # fixing module string
        if self.env.MODULE_STRING != '':
            match = re.search(r'START(.*?)EOF', self.env.MODULE_STRING, re.DOTALL)
            if match:
                extracted_string = match.group(1).strip()
                self.env.MODULE_STRING = extracted_string
            else:
                print('Can\'t find START and EOF to recognize the\
                              module string, using default...')
                self.env.MODULE_STRING = ''

        # checking exes
        if not os.path.isfile(self.env.LAMMPS_EXE):
            print('specified file {} not exist, using default...'.format(self.env.LAMMPS_EXE))
            self.env.LAMMPS_EXE = '/depot/bsavoie/apps/lammps/exe/lmp_mpi_180501'
        if not os.path.isfile(self.env.ORCA_EXE):
            print('specified file {} not exist, using default...'.format(self.env.ORCA_EXE))
            self.env.ORCA_EXE = '/depot/bsavoie/apps/orca_4_1_2/orca'

        if self.ff_fit.FUNCTIONAL != 'wB97X-D3':
            self.ff_fit.D3STR = 'D3BJ'
        else:
            self.ff_fit.D3STR = ''

        # don't try to overide the previous fitted db
        if os.path.isfile(PARAM_DB):
            raise IOConfigException('ERROR: {} already exist'.format(PARAM_DB))

        if not os.path.isfile(self.ff_fit.FF):
            print('WARNING: specified FF db {} does not exist'.format(self.ff_fit.FF))
            with open(PARAM_DB, 'w') as f:
                f.write('#emtpy db\n')
        else:
            shutil.copyfile(self.ff_fit.FF, PARAM_DB)

        self.ff_fit.FF = "{}/{}".format(self._running_dir, PARAM_DB)
    def ParseAndValidate(self, infile):
        self.Parse(infile)
        self.Validate()


class SubmitJobConfig:
    def __init__(self, **kwargs):
        # processes to use to single job
        self.PROCS = kwargs.get('procs', 4)
        # wall time for single job
        # for can accep min but need to be input as 30min (no space)
        self.WT = kwargs.get('wt', 4)
        # queue
        self.Q = kwargs.get('queue', 'bsavoie')
        # scheduler system
        self.SCHED = kwargs.get('sched', 'slurm-halstead')
        # ppn, related to machine
        self.PPN = kwargs.get('ppn', 128)
        # number of jobs to be bundled together
        self.SIZE = kwargs.get('size', 1)

    def FixType(self):
        self.PROCS = int(self.PROCS)
        self.PPN = int(self.PPN)
        self.SIZE = int(self.SIZE)

def main(argv):
    config = InputConfig()
    config.ParseAndValidate('test.config')
    print(config.param_ba.PROCS)
    config.WriteConfig('test.config2')


if __name__ == '__main__':
    main(sys.argv[1:])
import re
import os
import yaml as yml
from types import SimpleNamespace
# RE: PyYAML yaml.load(input) Deprecation (https://msg.pyyaml.org/load)
yml.warnings({'YAMLLoadWarning': False})

# set configuration filepaths
def _get_this_source_file_directory():
    def path(fullpath):
        return re.sub(r"(.*/).*$",r"\1",fullpath)
    return path(os.path.realpath(__file__))
config_main   = os.path.abspath(_get_this_source_file_directory()+"../config/config.yml")
build_context = os.path.abspath(_get_this_source_file_directory()+"../config/docker-compose.yml")

def cfg_get_main():
    return yml.load(open(config_main))

def cfg_get_build_context(module=None):
    context = yml.load(open(build_context))
    if not module is None:
        if module in context['services']:
            context = context['services'][module]
        else:
            context = None
    return context

# get survey beam pars from main configuration
def cfg_fetch_beam_pars(survey):
    if not survey is None:
        survey  = survey.lower()
        surveys = cfg_get_main()['surveys']
        if survey in surveys and 'beam' in surveys[survey]:
            beam = surveys[survey]['beam']
            try:
                return {'BMAJ': beam['bmaj'], 'BMIN': beam['bmin'], 'BPA': beam['bpa']} 
            except:
                pass
    return None

# get survey frequency from main configuration
def cfg_fetch_frequency(survey):
    if not survey is None:
        survey = survey.lower()
        surveys = cfg_get_main()['surveys']
        if survey in surveys and 'freq' in surveys[survey]:
            return surveys[survey]['freq']
    return None

#                  * * *   D E P R E C A T E D   * * *
# find if source finder/s require/s frequency as input
def is_frequency_required(module=None):
    modules = {list(m.keys())[0]:m[list(m.keys())[0]] for m in cfg_get_main()['modules']}
    if module is None:
        for key in modules:
            module = modules[key]
            if 'required' in module and not module['required'] is None and 'frequency' in module['required']:
                    return True
    elif module in modules:
        module = modules[module]
        if 'required' in module and 'frequency' in module['required']:
            return True
    return False

from .exceptions import HydraIOError
from .exceptions import HydraFreqError
from .exceptions import HydraConfigError
from .exceptions import print_warning
from .exceptions import print_ok
class Config:
    def __init__(self,fits_image_file=None):
        # make sure file exists before we start spawning containers
        self.image_file = fits_image_file
        if not fits_image_file is None and not os.path.isfile(self.image_file):
            msg = ""
            msg += f"> ERROR: File does not exist: {fits_image_file}"
            print_warning(msg)
            raise HydraIOError

        # loop over modules to build master config
        self.main    = cfg_get_main()
        self.modules = dict()
        for source_finder,configuration in {list(m.keys())[0]:m[list(m.keys())[0]] for m in self.main['modules']}.items():
            # debug
            #print(f"{source_finder}: {configuration}")

            # ok, let's build our module cfg
            module = dict()

            # assign name field -- must exist
            if 'name' in configuration:
                module['name'] = configuration['name']
            else:
                self.__print_error(f"Field {source_finder} 'name' missing: {config_main}")
                raise HydraConfigError

            # get container config file path and load
            if 'config' in configuration:
                cfg_container_file = os.path.abspath(f"{_get_this_source_file_directory()}/../config/{configuration['config']}")
                if not os.path.isfile(cfg_container_file):
                    self.__print_error(f"{module['name']} container config file not found: {cfg_container_file}")
                    raise HydraConfigError
                module['config'] = cfg_container_file
            else:
                self.__print_error(f"{module['name']} 'config' field missing: {config_main}")
                raise HydraConfigError
            cfg_container = yml.load(open(cfg_container_file))

            # get module credits
            module['credit'] = cfg_container['credit'] if 'credit' in cfg_container else None

            # get module input parameters
            if 'parameters' in cfg_container:
                parameters = cfg_container['parameters']
                module['variables']   = parameters['variables']   if 'variables' in parameters else None
                module['required']    = parameters['required']    if 'required'  in parameters else None
                #   * * *   T O   B E   D E P R E C A T E D   * * * 
                #module['suggested']   = parameters['suggested']   if 'suggested' in parameters else None
                module['constraints'] = parameters['constraints'] if 'constraints' in parameters else None
                module['optional'] = parameters['optional'] if 'optional' in parameters else None

            # debug
            #print(module)
            
            # get module output catalogue meta data
            if 'catalogue' in cfg_container:
                module['catalogue_meta'] = cfg_container['catalogue']
            else:
                self.__print_error(f"{module['name']} 'catalogue' definition block missing: {cfg_container_file}")
                raise HydraConfigError

            # push module info onto modules stack
            self.modules[source_finder] = module

            # debug
            #print("***")

        #print(f"Modules: {self.modules}")

    def __print_error(self,msg):
        print_warning(f"ERROR:{type(self).__name__}: {msg}")

    def __normalize(self,module,parameter):
        return "--"+module+"-"+re.sub(r"_","",parameter).lower()

    def get_required(self):
        required = list()
        for module,datum in self.modules.items():
            if 'required' in datum and not datum['required'] is None:
                for param in datum['required']:
                    #required.append(self.__normalize(module,list(param.keys())[0]))
                    required.append(list(param.keys())[0])
        return required

    def is_required(self,parameter):
        return parameter in self.get_required()

    def get_modules(self):
        return list(self.modules.keys())

    def get_variables(self):
        variables = dict()
        for module,datum in self.modules.items():
            if 'variables' in datum and not  datum['variables'] is None:
                variables[module] = datum['variables']
        return variables

    def get_parameters(self):
        def get_info(variable,param):
            info = param[variable]
            if  not 'default' in info:
                info['default'] = None
            if  not 'class' in info:
                info['class'] = None
            info['container_token'] = "--"+re.sub(r"_","-",variable)
            info['required'] = False
            info['optional'] = False
            return info
        variables = list()
        for module,datum in self.modules.items():
            if 'variables' in datum and not datum['variables'] is None:
                for param in datum['variables']:
                    variable = list(param.keys())[0]
                    #info = param[variable]
                    #if  not 'default' in info:
                    #    info['default'] = None
                    #info['container_token'] = "--"+re.sub(r"_","-",variable)
                    #info['required'] = False
                    #info['optional'] = False
                    #info['class'] = None
                    info = get_info(variable,param)
                    variables.append({self.__normalize(module,variable): info})
            if 'required' in datum and not datum['required'] is None:
                for param in datum['required']:
                    variable = list(param.keys())[0]
                    #info = param[variable]
                    #info['default'] = None
                    #info['required'] = True
                    #info['optional'] = False
                    #info['container_token'] = "--"+re.sub(r"_","-",variable)
                    info = get_info(variable,param)
                    info['required'] = True
                    variables.append({self.__normalize(module,variable): info})
            if 'optional' in datum and not datum['optional'] is None:
                for param in datum['optional']:
                    variable = list(param.keys())[0]
                    #info = param[variable]
                    #if  not 'default' in info:
                    #    info['default'] = None
                    #if  not 'class' in info:
                    #    info['class'] = None
                    #info['container_token'] = "--"+re.sub(r"_","-",variable)
                    #info['required'] = False
                    #info['optional'] = True
                    info = get_info(variable,param)
                    info['optional'] = True
                    variables.append({self.__normalize(module,variable): info})
        return variables

    def get_constraints(self):
        constraints = dict()
        for module,datum in self.modules.items():
            info = dict()
            if 'constraints' in datum and not datum['constraints'] is None:
                if 'range' in datum['constraints'] and not datum['constraints']['range'] is None:
                    info['parameters'] = datum['constraints']['range']
                    defaults = dict()
                    if 'variables' in datum and not datum['variables'] is None:
                        for param in info['parameters']:
                            for item in datum['variables']:
                                variable = list(item.keys())[0] 
                                if param == variable and 'default' in item[variable] and not item[variable]['default'] is None:
                                    defaults[variable] = item[variable]['default']
                                    break
                    info['defaults'] = defaults
                else:
                    info['parameters'] = dict()
                if 'directive' in datum['constraints'] and not datum['constraints']['directive'] is None:
                    info['directive'] = datum['constraints']['directive']
                else:
                    info['directive'] = True
            else:
                info['parameters'] = None
                info['directive']  = None
            constraints[module] = info
        return constraints


    def has_rms_box_pars(self,module):
        if module in self.modules and not self.modules[module]['optional'] is None:
            rms_box_pars = list()
            for datum in self.modules[module]['optional']:
                param = list(datum.keys())[0]
                if 'class' in datum[param] and 'rms_box' in datum[param]['class']:
                    rms_box_pars.append(param)
            if 'box_size' in rms_box_pars and 'step_size' in rms_box_pars:
                return True
        return False


    def get_catalogue_meta(self,module=None):
        if module is None:
            meta = dict()
            for module in self.modules:
                meta[module] = self.modules[module]['catalogue_meta']
            return meta
        return self.modules[module]['catalogue_meta']

    def get_catalogue_modules(self):
        return self.main['modules']

    def get_catalogue_common(self):
        return self.main['common']

    def get_graphics_context(self):
        return self.main['graphics']






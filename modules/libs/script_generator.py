import re
import sys
sys.path.insert(0,"../")
from jinja2 import Template
from libs.config_tools import Config
from datetime import date

if __name__ == "__main__":
    cfg = Config()
    variables = cfg.get_parameters()
    #print(variables)

    with open('templates/cerberus.py-jn2') as fd:
        template = Template(fd.read())
    with open('cerberus.py','w') as fd:
        def mk_var(param):
            return re.sub(r"-","_",re.sub(r"^-+","",param))
        env = list()
        modules = list()
        for variable in variables:
            flag = list(variable.keys())[0]
            raw_flag = re.sub(r"^-+.*?-",r"--",flag)
            data = variable[flag]
            module = re.sub(r"^-+(.*?)-.*$",r"\1",flag)
            # debug
            #var=re.sub(r"-","_",re.sub(r"^-+","",flag))
            #rvar=re.sub(r"^-+","",raw_flag)
            #print(f"{module}: flag={flag} (var={var}), raw_flag={raw_flag} (rvar={rvar})")
            env.append({
                'module': module,
                'flag': flag,
                'rflag': raw_flag,
                'token': data['container_token'],
                'var': re.sub(r"-","_",re.sub(r"^-+","",flag)),
                'rvar': re.sub(r"^-+","",raw_flag),
                'type': data['type'],
                'default': data['default'],
                'required': data['required'],
                'optional': data['optional'],
                'class': data['class'],
                'help': data['definition'],
            })
            modules.append(module)
        fd.write(template.render(modules=sorted(set(modules)),variables=env,date=date.today()))

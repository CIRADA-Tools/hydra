class HydraIOError(Exception):
    pass

class HydraFreqError(Exception):
    pass

class HydraConfigError(Exception):
    pass

class HydraContainerError(Exception):
    pass

class HydraRenderingError(Exception):
    pass

# notes: https://stackoverflow.com/questions/287871/how-to-print-colored-text-in-python
class bcolors:
    OKBLUE    = '\033[94m'
    OKGREEN   = '\033[92m'
    WARNING   = '\033[91m'
    ENDC      = '\033[0m'
    BOLD      = '\033[1m'
    UNDERLINE = '\033[4m'

def print_warning(msg):
    print(f"{bcolors.WARNING}{bcolors.BOLD}{msg}{bcolors.ENDC}")

def print_ok(msg):
    print(f"{bcolors.OKBLUE}{bcolors.BOLD}{msg}{bcolors.ENDC}")


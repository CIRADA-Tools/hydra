import pyximport
pyximport.install(language_level=3)
from .renderizer import renderer

def render(datum):
    renderer(datum)

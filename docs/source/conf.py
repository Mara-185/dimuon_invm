# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
package_name = 'dimuon_invm'
package_name2 = 'utils'
package_name3 = 'Z_asymmetry'
package_name4 = 'feynman_Z'
package_root = os.path.abspath('../..')
#sys.path.insert(0, package_root)
sys.path.insert(0, os.path.join(package_root, package_name))


package_root2 = os.path.abspath('../../Utils/')
sys.path.insert(0, package_root2)
sys.path.insert(0, os.path.join(package_root2, package_name2))

#package_root3 = os.path.abspath('../../Z_asymmetry/')
sys.path.insert(0, os.path.join(package_root, package_name3))
sys.path.insert(0, os.path.join(package_root2, package_name4))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'dimuon_invm'
copyright = "2022, Mara Stefania Calo'"
author = "Mara Stefania Calo'"
release = '0.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
     'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',]

templates_path = ['_templates']
exclude_patterns = [u'_build', 'Thumbs.db', '.DS_Store']

source_suffix = '.rst'
master_doc = 'index'
pygments_style = 'sphinx'

autodoc_mock_imports = ['ROOT']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

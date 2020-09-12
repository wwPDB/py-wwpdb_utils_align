# File: setup.py
# Date: 14-Oct-2018 jdw
#
# Update:
#
import glob
import os
import platform
import re
import subprocess
import sys
from distutils.version import LooseVersion

from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir="", sources=None):
        sources = sources if sources else []
        Extension.__init__(self, name, sources=sources)
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " + ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmakeVersion = LooseVersion(re.search(r"version\s*([\d.]+)", out.decode()).group(1))
            if cmakeVersion < "3.1.0":
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):

        #
        debug = True
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmakeArgs = ["-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + extdir, "-DPYTHON_EXECUTABLE=" + sys.executable]

        # We need to help cmake find the correct python for this virtual env -
        # ---
        libPath = None
        lsp = os.path.join(sys.exec_prefix, "lib", "libpython") + "*"
        lpL = glob.glob(lsp)
        if lpL:
            libPath = lpL[0]
        elif hasattr(sys, "base_exec_prefix"):
            lsp = os.path.join(sys.base_exec_prefix, "lib", "libpython") + "*"  # pylint: disable=no-member
            lpL = glob.glob(lsp)
            if lpL:
                libPath = lpL[0]
        if libPath:
            cmakeArgs += ["-DPYTHON_LIBRARY=" + libPath]
        else:
            print("------ WARNING could not locate python library")
        # ---
        inclPath = None
        isp = os.path.join(sys.exec_prefix, "include", "python") + "%s.%s" % (sys.version_info.major, sys.version_info.minor) + "*"
        ipL = glob.glob(isp)
        if ipL:
            inclPath = ipL[0]
        elif hasattr(sys, "base_exec_prefix"):
            isp = os.path.join(sys.base_exec_prefix, "include", "python") + "%s.%s" % (sys.version_info.major, sys.version_info.minor) + "*"  # pylint: disable=no-member
            ipL = glob.glob(isp)
            if ipL:
                inclPath = ipL[0]
        if inclPath:
            cmakeArgs += ["-DPYTHON_INCLUDE_DIR=" + inclPath]
        else:
            print("------ WARNING could not locate python include files")
        # ---
        cfg = "Debug" if self.debug else "Release"
        buildArgs = ["--config", cfg]

        if platform.system() == "Windows":
            cmakeArgs += ["-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(cfg.upper(), extdir)]
            if sys.maxsize > 2 ** 32:
                cmakeArgs += ["-A", "x64"]
            buildArgs += ["--", "/m"]
        else:
            cmakeArgs += ["-DCMAKE_BUILD_TYPE=" + cfg]
            buildArgs += ["--", "-j2"]

        env = os.environ.copy()
        env["CXXFLAGS"] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get("CXXFLAGS", ""), self.distribution.get_version())
        env["RUN_FROM_DISUTILS"] = "yes"
        #
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        #
        if debug:
            print("------------- setup.py -----------------")
            print("Extension source path ", ext.sourcedir)
            print("CMAKE_ARGS ", cmakeArgs)
            print("self.build_temp ", self.build_temp)
            print("extdir", extdir)
            print("ext.name", ext.name)
            print("sys.executable", sys.executable)
            print("sys.exec_prefix", sys.exec_prefix)
            print("CXXFLAGS ", env["CXXFLAGS"])

        #
        subprocess.check_call(["cmake", ext.sourcedir] + cmakeArgs, cwd=self.build_temp, env=env)
        subprocess.check_call(["cmake", "--build", "."] + buildArgs, cwd=self.build_temp)


packages = []
thisPackage = "wwpdb.utils.align"
requires = []


with open("wwpdb/utils/align/__init__.py", "r") as fd:
    version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]', fd.read(), re.MULTILINE).group(1)

if not version:
    raise RuntimeError("Cannot find version information")

setup(
    name=thisPackage,
    version=version,
    description="Alignment Library and Tools",
    long_description="See:  README.md",
    author="John Westbrook",
    author_email="john.westbrook@rcsb.org",
    url="http://github.com/wwpdb/py-wwpdb_utils_align",
    #
    license="Apache 2.0",
    classifiers=[
        "Development Status :: 3 - Alpha",
        # 'Development Status :: 5 - Production/Stable',
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
    # entry_points={
    #    'console_scripts': [
    #        'onedep_validate_cli=onedep.cli.validate_cli:run',
    #    ]
    # },
    #
    install_requires=["future", "six"],
    packages=find_packages(exclude=["wwpdb.utils.tests-align", "tests.*"]),
    package_data={
        # If any package contains *.md or *.rst ...  files, include them:
        "": ["*.md", "*.rst", "*.txt", "*.h", "*.C", ".c", "*.cpp"],
    },
    #
    #
    test_suite="wwpdb.utils.tests-align",
    tests_require=["tox"],
    #
    # Not configured ...
    extras_require={
        "dev": ["check-manifest"],
        "test": ["coverage"],
    },
    # Added for
    command_options={"build_sphinx": {"project": ("setup.py", thisPackage), "version": ("setup.py", version), "release": ("setup.py", version)}},
    ext_modules=[CMakeExtension("wwpdb.utils.align.alignlib")],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
)

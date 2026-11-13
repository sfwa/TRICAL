import os
import shutil
import subprocess

from setuptools import setup
from setuptools.command.build_py import build_py

ROOT = os.path.dirname(os.path.abspath(__file__))
LIBRARY_NAMES = ("libTRICAL.so", "libTRICAL.dylib", "TRICAL.dll",
                 "libTRICAL.dll")


class build_trical(build_py):
    description = """Build the TRICAL shared library"""

    def run(self):
        # Build the shared library out of tree, then put it in the package
        build_dir = os.path.join(ROOT, "build", "trical")
        os.makedirs(build_dir, exist_ok=True)

        subprocess.check_call(["cmake", "-DCMAKE_BUILD_TYPE=Release", ROOT],
                              cwd=build_dir)
        subprocess.check_call(["cmake", "--build", ".", "--target", "TRICAL",
                               "--config", "Release"], cwd=build_dir)

        package_dir = os.path.join(ROOT, "python", "TRICAL")
        for directory, _, files in os.walk(build_dir):
            for name in files:
                if name in LIBRARY_NAMES:
                    shutil.copy(os.path.join(directory, name), package_dir)

        build_py.run(self)


with open(os.path.join(ROOT, "README.md")) as f:
    long_description = f.read()

setup(
    name="trical",
    url="https://github.com/sfwa/TRICAL",
    author="Ben Dyer",
    author_email="",
    version="2.0.0",
    description="""Real-time scale and bias calibration algorithm for \
tri-axial field sensors (e.g. magnetometers).""",
    long_description=long_description,
    long_description_content_type="text/markdown",
    package_dir={"": "python"},
    packages=["TRICAL"],
    package_data={
        "TRICAL": list(LIBRARY_NAMES) + ["viz.html"]
    },
    license="MIT License",
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Topic :: Software Development :: Libraries"
    ],
    python_requires=">=3.6",
    cmdclass={"build_py": build_trical}
)

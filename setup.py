from setuptools import setup, find_packages

setup(
    name='antigen-antibody-residue-hotspot', 
    version='0.1.0',
    author='Your Name',
    author_email='your.email@example.com',
    description='A tool for antigen-antibody residue hotspot analysis using Ridge regression.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/<your-username>/antigen-antibody-residue-hotspot',
    packages=find_packages(),
    license='MIT',
    install_requires=[
        'numpy',
        'pandas',
        'scikit-learn',
        'matplotlib',
        'seaborn'
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)

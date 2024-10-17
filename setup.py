import setuptools

version = [l.strip() for l in open("malaria_profiler/__init__.py") if "version" in l][0].split('"')[1]

setuptools.setup(

	name="malaria_profiler",

	version=version,
	packages=["malaria_profiler"],
	license="GPLv3",
	long_description="Malaria-Profiler command line tool",
	scripts= [
		'scripts/malaria-profiler',
        'scripts/mp-analyse-nanopore-vivax-amplicons.py',
        'scripts/mp-collate-results.py'
		],
)

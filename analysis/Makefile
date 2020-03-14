PY = $$PYTHONEXE
DATA = ../data

all: test

test:
#	nosetests -v
	$(PY) wrightfisher_test.py

prof:
	$(PY) -m cProfile -s time timewf.py

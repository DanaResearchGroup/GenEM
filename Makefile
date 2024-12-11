test:
	python -m unittest discover

run-iso:
	python -m src.isostere_optimization

run-func:
	python -m src.functional_group_optimization

run:
	python -m src.main
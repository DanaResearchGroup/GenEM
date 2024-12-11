test:
	python -m unittest discover

run-back:
	python -m src.backbone_optimization

run-func:
	python -m src.functional_group_optimization

run-atom:
	python -m src.atom_substitution_optimization

run:
	python -m src.main
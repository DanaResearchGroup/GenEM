test:
	python -m unittest discover

run-back:
	python -m src.backbone_optimization

run-func:
	python -m src.functional_group_optimization

run:
	python -m src.main
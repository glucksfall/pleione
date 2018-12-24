push-github:
	git add .
	git commit -a
	git push --repo=git@github.com:glucksfall/pleione.git

publish-pypi:
	versioneer install -f
	python3 setup.py sdist bdist_wheel
	twine upload dist/*

.PHONY: all figures pdf markdown test clean

all: figures pdf test

figures:
	python3 scripts/generate_all.py

pdf: figures
	latexmk -xelatex -interaction=nonstopmode -halt-on-error -cd notes/brownian_motion_notes.tex

markdown:
	pandoc notes/brownian_motion_notes.tex --from=latex --to=gfm --wrap=none --output=notes/brownian_motion_notes.md
	python3 scripts/normalize_markdown.py

test:
	python3 -m unittest discover -s tests -v

clean:
	latexmk -c -cd notes/brownian_motion_notes.tex

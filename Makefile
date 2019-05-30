#
JUPYTER = jupyter

vpath %.ipynb notebooks

.SUFFIXES: .md .ipynb

$(OUTDIR)/%.md: %.ipynb
	PYTHONPATH=ipython/convert/ $(JUPYTER) nbconvert --to markdown --output-dir $(OUTDIR) --config ipython/convert/ipython_nbconvert_config.py $<

.PHONY: clean

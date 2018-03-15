

archive_examples_data.zip : examples/data/
	zip -r $@ $<

.PHONY: unpack
unpack:
	unzip archive_examples_data.zip

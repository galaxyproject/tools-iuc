planemo tool_init --force \
				   --id 'igblast-parser' \
				   --name 'convert igblast output into simple csv table' \
				   --requirement igblast-parser@0.0.3\
				   --example_command 'igblast-parser --in <igblast.output> --out <parser_output>' \
				   --example_input test-data/igblast_output.txt \
				   --example_output test-data/parser_output.csv \
				   --test_case \
				   --cite_url 'https://github.com/aerijman/igblast-parser' \
				   --help_from_command 'igblast-parser --help'

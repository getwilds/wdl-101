
# wdl-101
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)

This is the repository for the WDL and JSON files accompanying the [WDL 101 course](https://hutchdatascience.org/WDL_Workflows_Guide/index.html).
You will find here the following:
1. **mutation_calling.wdl**: this is a working WDL workflow
2. **mutation_calling_input_fh.json**: this JSON is meant for use by the Fred Hutch community. It points to file paths that are accessible only with Fred Hutch credentials
3. **MOLM13_combined_final.fastq, HCC4006_final.fastq, CALU1_combined_final.fastq**: These are the fastq files that can be used as inputs in the WDL. If you do not have Fred Hutch credentials and want to use these files to test this or any other workflow download a copy of this locally and make sure to edit the input.JSON file with the full path to where these files are stored.
4. **mutation_calling_input.json**: This would be the input JSON you would want to enter the full file paths to the local copy of the sample files referenced above. 

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on our [issue tracker](https://github.com/getwilds/fastq-to-cram/issues).

## Contributing

If you would like to contribute to this WILDS WDL workflow, see our [contribution guidelines](.github/CONTRIBUTING.md) as well out our [WILDS Contributor Guide](https://getwilds.org/guide/) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.

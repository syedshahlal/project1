# Mini-Project 1: Parsing Variant Call Format File

## Introduction

This is the first of two mini-projects in the masters course. In this project, we will work with a variant call format file containing one thousand variants. Your task is to parse each line of the file, transform it into dictionaries, and generate a list of dictionaries as the final output.

### Data Source

The variant call format file can be found [here](https://en.wikipedia.org/wiki/Variant_Call_Format). Each line in the file represents a variant. The file contains header lines, which start with one hash symbol (#), followed by variant lines.

## Task Overview

Your goal is to create Python functions to perform the following tasks:

1. Parse the variant file and extract specific predictor fields.
2. Convert text descriptions of predictors to integer values and calculate a sum.
3. Process a gzipped variant file using the provided code snippet.
4. Select variants with non-zero sum predictor values and save the output.

### Dictionary Structure

The output of this project will be a list of dictionaries, where each dictionary will have the following structure:

```python
{
    "ALT": "G",
    "CHROM": "4",
    "FILTER": "PASS",
    "ID": ".",
    "INFO": {
        # ... Predictor information ...
    },
    "POS": 123416186,
    "QUAL" : 23.25,
    "REF": "A",
    "SAMPLE": {
        "XG102": {
            # ... Sample information ...
        }
    }
}
```

## Task Breakdown

### Parsing and Extracting Fields

1. Read `mini_project_data.json`.
2. Select variants with specific predictor values.
3. Extract fields like `CHROM`, `POS`, `REF`, and `ALT`.

### Converting Predictor Texts to Integers

1. Convert predictor text descriptions to integer values.
2. Sum up the integer values for each variant.

### Processing Gzipped Files

1. Use the provided code snippet to read a gzipped file.
2. Repeat the parsing and field extraction process for the gzipped file.

### Selecting Variants with Non-Zero Predictors

1. Open `mini_project1_gzip.json`.
2. Select variants with non-zero sum predictor values.

## Implementation Details

The provided description outlines the high-level steps for each function. Follow the instructions to complete each function and generate the required outputs.

## Additional Resources

For additional resources, refer to the provided links and documentation:

- [Variant Call Format](https://en.wikipedia.org/wiki/Variant_Call_Format)
- [Predictor Fields](https://brb.nci.nih.gov/seqtools/colexpanno.html#dbnsfp)

Happy coding!

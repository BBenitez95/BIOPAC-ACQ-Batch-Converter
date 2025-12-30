# BIOPAC ACQ Batch Converter

A Shiny app for batch converting BIOPAC AcqKnowledge (.acq) files to text (.txt) format.

![R](https://img.shields.io/badge/R-%3E%3D4.0-blue)
![License](https://img.shields.io/badge/license-MIT-green)

## Features

- Batch convert multiple .acq files or entire folders
- Preview signal data with interactive multi-channel plots
- Export to tab-delimited (.txt) or CSV (.csv) format
- Optional metadata headers
- Supports AcqKnowledge versions 3.x through 5.x (uncompressed)

## Installation

```r
install.packages("pacman")
pacman::p_load(shiny, shinyjs, shinyFiles, plotly)
```

## Usage

```r
shiny::runApp("acq_batch_converter.R")
```

1. Select input files or folder containing .acq files
2. Click **Convert** to process files
3. Review signals and data in the preview panels
4. Select output folder (defaults to input location)
5. Click **Save** to export all converted files

## Output Format

**With metadata:**
```
# Source: recording.acq
# Version: 4.4.0
# Sample Rate: 2000 Hz
# Channels: 3
# Samples: 10000
#
# Ch1: EMG (mV)
# Ch2: Force (N)
# Ch3: Trigger (V)
#
time_sec	EMG	Force	Trigger
0.0000	0.001234	12.345678	0.000000
0.0005	0.001456	12.346789	0.000000
```

**Without metadata:**
```
time_sec	EMG	Force	Trigger
0.0000	0.001234	12.345678	0.000000
0.0005	0.001456	12.346789	0.000000
```

## Limitations

- Compressed .acq files are not supported (resave as uncompressed in AcqKnowledge)
- Versions prior to 3.0 are untested

## Acknowledgments

ACQ binary parser ported from [bioread](https://github.com/uwmadison-chm/bioread) by Nate Vack, UW-Madison.

## License

MIT

# BiteID - Individual Identification via Dental Models

BiteID is a 3D Slicer module that enables individual identification based on 3D dental models of the mandible or maxilla. It works by automatically registering two dental scans, establishing point correspondences between them, and computing metrics of overall similarity for individual identification purposes.

## Features

- Perform dental scan registration and similarity calculation
- Utilizes open3d for point cloud processing
- Configurable parameters for advanced processing options
- Batch processing of multiple models against a model database 
- Live preview of ID matching results

## Requirements

- 3D Slicer version 5.2.2 or later
- open3d version 0.17 (the module will prompt you for installation if it's not found)

## Installation

1. **Clone the repository**: The first step would be to download the BiteID repository to your local machine. 

2. **Launch 3D Slicer and open the Extension Wizard**: Start 3D Slicer on your computer. In the menu bar, go to "Find Module" (looking glass) -> "Extension Wizard". 

3. **Add the cloned repository**: In the Extension Wizard, there should be an option to add a new extension from a local directory. Select this, then navigate to the directory where you cloned the BiteID repository.

4. **Configure and generate the extension**: The Extension Wizard will examine the contents of the repository and should provide you with a list of options to configure the build. 

Please note that the module is currently available under the 'Examples'. This will be changed in the future.

## Usage

- Load 3D dental models into 3D Slicer
- Open the BiteID module in 3D Slicer
- Select the source and target models for comparison
- Run the BiteID module
- Adjust settings for better results if necessary


## License

This project is licensed under the BSD-2 License.

## Acknowledgements


## Contributions

Contributions to this project are welcome. Please open a pull request with your changes, or open an issue to discuss what you would like to change.

## Questions

If you have any questions or issues, please open an [issue](https://github.com/agporto/BiteID/issues).

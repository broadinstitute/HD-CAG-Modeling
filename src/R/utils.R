
# This must be pre-loaded (i.e. not deferred until CAGModeling$init)
if (!exists("R6Class")) {
    library(R6, quiet=TRUE, warn=FALSE)
}

CAGModeling <- R6Class("CAGModeling", list())
CAGModeling$initialized = FALSE
CAGModeling$suppress.messages = TRUE

CAGModeling$loadLibrary = function(pkg) {
    if (CAGModeling$suppress.messages) {
        suppressWarnings(suppressMessages(library(pkg, character.only=TRUE, quiet=TRUE, warn=FALSE)))
    } else {
        library(pkg, character.only=TRUE)
    }
}

CAGModeling$init = function(load.philentropy=FALSE, load.optimization=FALSE) {

    if (!CAGModeling$initialized) {
        CAGModeling$loadLibrary("R6")
        CAGModeling$loadLibrary("data.table")
        CAGModeling$loadLibrary("tidyverse")
        CAGModeling$loadLibrary("dplyr")
        CAGModeling$loadLibrary("pracma")
        CAGModeling$loadLibrary("reshape")
        CAGModeling$loadLibrary("cowplot")
        CAGModeling$loadLibrary("RColorBrewer")
        CAGModeling$loadLibrary("ggfortify")
        CAGModeling$loadLibrary("expm")
        CAGModeling$initialized = TRUE
    }

    # Optional libraries (not currently used)
    if (load.philentropy) {
        CAGModeling$loadLibrary("philentropy")
    }
    if (load.optimization) {
        CAGModeling$loadLibrary("optimization")
    }
    return(invisible(NULL))
}

CAGModeling$loadModels = function(inputFile) {
    return(CAGBase$load(inputFile))
}

CAGModeling$saveModels = function(modelList, outputFile) {
    outputData = CAGBase$encode_value(modelList)
    saveRDS(outputData, outputFile)
}

CAGModeling$mergeModelFiles = function(inputFileList, outputFile) {
    outputDataList = NULL
    for (inputFile in inputFileList) {
        modelData = readRDS(inputFile)
        if (!is.null(modelData)) {
            if (CAGBase$is_object_state(modelData)) {
                outputDataList = append(outputDataList, list(modelData))
            } else {
                outputDataList = append(outputDataList, modelData)
            }
        }
    }
    saveRDS(outputDataList, outputFile)
}

CAGModeling$readDonorInfo = function(fileName) {
    result = NULL
    fileData = fread(fileName, header=T, sep="\t")
    donors = fileData$Donor
    if (length(unique(donors)) != length(donors)) {
        warn(sprintf("Duplicate donor names in input file: %s", fileName))
        return(NULL)
    }
    for (donor in donors) {
        donorData = fileData[fileData$Donor == donor]
        result = c(result, CAGDonor$new(donor, age=donorData$Age, inherited_cag=donorData$CAG))
    }
    names(result) = donors
    return(result)
}

CAGModeling$readCellLoss = function(fileName, dataSetList) {
    fileData = fread(fileName, header=T, sep="\t")
    for (dataSet in dataSetList) {
        donor = dataSet$donor$id
        region = dataSet$region
        celltype = dataSet$celltype
        loss = fileData[fileData$DONOR == donor & fileData$REGION == region & fileData$CELLTYPE == celltype]$CELL_LOSS_FRACTION
        if (!is.null(loss) && length(loss) > 0) {
            dataSet$loss_fraction = loss
        }
    }
    return(dataSetList)
}

CAGModeling$readCellFile = function(fileName, donorMap, donors=NULL, regions=NULL, celltypes=NULL, minCAG=35) {
    result = NULL
    fileData = fread(fileName, header=T, sep="\t")
    cagData = fileData[fileData$REPLENGTH >= minCAG]
    donorList = sort(unique(cagData$DONOR))
    if (!is.null(donors)) {
        donorList = intersect(donorList, donors)
    }
    donorList = intersect(donorList, names(donorMap))
    for (donor in donorList) {
        donorInfo = donorMap[[donor]]
        donorData = cagData[cagData$DONOR == donor]
        regionList = sort(unique(donorData$REGION))
        if (!is.null(regions)) {
            regionList = intersect(regionList, regions)
        }
        for (region in regionList) {
            regionData = donorData[donorData$REGION == region]
            cellTypeList = sort(unique(regionData$CELLTYPE))
            if (!is.null(celltypes)) {
                cellTypeList = intersect(cellTypeList, celltypes)
            }
            for (cellType in cellTypeList) {
                cellData = regionData[regionData$CELLTYPE == cellType]
                cags = as.integer(round(cellData$REPLENGTH))
                cags = sort(cags)
                if (length(cags) > 0) {
                    dataSet = CAGDataSet$new(donor=donorInfo, region=region, celltype=cellType, observed=cags)
                    result = c(result, dataSet)
                }
            }
        }
    }
    return(result)
}

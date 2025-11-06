## Load packages
library(binaryLogic)
library(stringr)
library(flowCore)

# library(sodium)


f.createFolder <- function(folderPath) {
  if (dir.exists(folderPath) == F) {
    dir.create(folderPath, showWarnings = F, recursive = T)
  }
}


f_open_folder <- function(dir = getwd()) {
  if (.Platform["OS.type"] == "windows") {
    shell.exec(dir)
  } else {
    system(paste(Sys.getenv("R_BROWSER"), dir))
  }
}

f.results.table <- function(folder.path) {
  filenames <- list.files(folder.path, pattern = ".html", full.names = T)
  filenames <- gsub(".html", "", filenames)
  
  file_urls <- paste0("file://", filenames)
  
  
print(      paste0(
  " <a href = '", filenames[1], ".html' style='color:black;font-weight:bold' target = '_blank' >",
  filenames, "</a>"
))
  data.frame(
    `Patient Reports` =
      paste0(
        " <a href = '", file_urls, ".html' style='color:black;font-weight:bold' target = '_blank' >",
        filenames, "</a>"
      )
  )
}


#' Functions to get the part of string after string split
#' names: a vector contains strings
#' i the part of the strings you want to extract
StrExtract <- function(strings, split = "_", i, supresswarnings = T) {
  y <- strsplit(strings, split)
  len <- unlist(lapply(y, length))
  if (length(which(len < i)) > 0 & supresswarnings == F) {
    warning("Strings cannot be extracted in some items")
  } else {}

  y <- lapply(y, function(x) {
    if (i == 0) {
      return(x[length(x)])
    } else if (length(x) < i) {
      ori <- paste(x, collapse = split)
      return(ori)
    } else {
      return(x[i])
    }
  })
  y <- unlist(y)
  return(y)
}

f.createFolder <- function(folderPath) {
  if (dir.exists(folderPath) == F) {
    dir.create(folderPath, showWarnings = F, recursive = T)
  }
}



#' Extract Channel and Marker Information from a flowFrame
#'
#' This function takes a flowFrame object and extracts detailed information about
#' its channels, including channel names, descriptions, fluorochromes, markers,
#' and channel types (scatter, time, fluorescence).
#'
#' @param ff A flowFrame object containing flow cytometry data
#' @return A data frame with channel information including:
#'   - name: Channel name
#'   - desc: Channel description
#'   - fluorochome: Extracted fluorochrome (from description after second "_")
#'   - marker: Extracted marker (from description before first "_")
#'   - dmf: Display-friendly name (combination of channel name and description)
#'   - type: Channel classification (scatter_forward, scatter_side, time, fluorescence, unknown)
#' @examples
#' \dontrun{
#'   library(flowCore)
#'   # Load a FCS file
#'   ff <- read.FCS("myfile.fcs")
#'   # Get channel information
#'   channel_info <- f.get_channels_and_markers(ff)
#' }
f.get_channels_and_markers <- function(ff) {
  # Extract parameter information from the flowFrame
  params <- parameters(ff)
  
  # Create a dataframe with the first two columns (name and description)
  result <- as.data.frame(params@data)[,1:2]
  
  # Extract fluorochrome information from the description
  # Assumes description format has fluorochrome after second "_"
  result$fluorochome <- StrExtract(result$desc, "_", 2)
  
  # Extract marker information from the description
  # Assumes description format has marker name before first "_"
  result$marker <-  StrExtract(result$desc, "_", 1)
  
  # Initialize display-friendly name column with just the channel name
  result$dmf <-  result$name 
  
  # Identify which channels have descriptions (stained channels)
  stained <- which(!is.na(result$desc))
  
  # For stained channels, create a display-friendly name combining channel name and description
  result$dmf[stained] <- paste0(result$name[stained], " ", result$desc[stained])
  
  # Classify channels based on naming patterns
  result$type <- "unknown"  # Default type
  
  # Identify forward scatter channels (FSC)
  result$type[grep("^FSC", result$name, ignore.case = TRUE)] <- "scatter_forward"
  
  # Identify side scatter channels (SSC)
  result$type[grep("^SSC", result$name, ignore.case = TRUE)] <- "scatter_side"
  
  # Identify time channels
  result$type[grep("^Time", result$name, ignore.case = TRUE)] <- "time"
  
  # Identify fluorescence channels - using pattern of three numbers followed by "-A"
  # This typically matches fluorescence detector channels like "123-A"
  result$type[grep("[0-9]{3}-A$", result$name, ignore.case = TRUE)] <- "fluorescence"
  
  # Return the complete channel information dataframe
  return(result)
}



f.pop_MFI <- function(ff) {
  
  info <- f.get_channels_and_markers(ff)
  chnls <-which(info$type == "fluorescence")
  res <- robustbase::colMedians(exprs(ff)[, chnls])
  names(res) <- info$dmf[chnls]
  
  res <- c(res, count = nrow(ff))
  
  
  
  
}



#' Order Control Files Based on Prefix and Numeric Value
#'
#' This function takes a vector of filenames and orders them based on a specific pattern:
#' a prefix (alphabetic characters) followed by a numeric value. The prefixes are sorted
#' according to a custom order (`UV`, `V`, `B`, `YG`, `R`), and the numeric values within
#' each prefix group are sorted from low to high.
#'
#' @param filenames A character vector of filenames containing the pattern:
#'   - A prefix (1-5 alphabetic characters, e.g., "UV", "B", "YG")
#'   - Followed by 3 digits (e.g., "537", "809", "730")
#'   Example: `c("B537", "UV809", "YG730")`
#'
#' @return An integer vector of indices that would sort the input `filenames` according to:
#'   - The custom prefix order (`UV`, `V`, `B`, `YG`, `R`)
#'   - The numeric value within each prefix group (low to high)
#'
#' @examples
#' \dontrun{
#'   # Example input
#'   filenames <- c("B537", "V470", "B710", "YG730", "UV809", "YG660", "B602")
#'
#'   # Get the order of filenames
#'   order_indices <- f.order_control_files(filenames)
#'
#'   # Apply the order to the filenames
#'   sorted_filenames <- filenames[order_indices]
#'   print(sorted_filenames)
#' }
#'
#' @importFrom stringr str_extract
f.order_control_files <- function(filenames) {
  # Extract the pattern: a space, followed by 1-5 alphabetic characters, then 3 digits, and another space
  pmt <- str_extract(filenames, " [A-Za-z]{1,5}\\d{3} ")
  pmt <- gsub(" ", "", pmt)
  # Extract the alphabetic prefix from the extracted pattern
  prefixes <- str_extract(pmt, "^[A-Za-z]+")
  
  # Extract the numeric part from the extracted pattern and convert it to a numeric value
  numbers <- as.numeric(sub("[A-Za-z]+", "", pmt))
  
  # Define the custom order for prefixes
  custom_order <- c("UV", "V", "B", "YG", "R")
  
  # Convert prefixes to a factor with the custom order
  prefixes <- factor(prefixes, levels = custom_order)
  
  # Sort the data: first by prefix order, then by numeric value
  order(prefixes, numbers)
}




#' Calculate Compensation Matrix
#'
#' This function calculates the compensation matrix for a set of flow cytometry data. It uses positive and negative beads to determine the compensation factors for each channel.
#'
#' @param gs A list of flow cytometry data frames.
#' @param beads_type The type of beads to use for the compensation calculation. Default is "beads".
#' @param spectral Logical. If TRUE, the function will return all fluoresent coloumns. Default is FALSE.
#'
#' @return A list with two elements:
#' - `matrix`: A compensation matrix as a data frame.
#' - `event_counts`: A named vector of event counts for each channel.
#'
#'
#' @examples
#' # Example usage
#' comp_matrix <- Calculate_Compensation_Matrix(my_flow_data, beads_type = "my_beads")
#'
Calculate_Compensation_Matrix <- function(gs, beads_type= "beads", spectral = F ) {
  
  # Get the flowframe
  pos_beads <- paste0("/", beads_type, "/pos")
  pos_beads_mfi <- lapply(gs, function(g) {
    ff <- gs_pop_get_data(g, pos_beads)
    f.pop_MFI(ff[[1]])
  })
  
  event_conts <- unlist(lapply(pos_beads_mfi, function(x) {
    x["count"]
  }))
  names(event_conts) <- StrExtract( names(event_conts), " ", 3)
  event_conts <- event_conts[-length(event_conts)]
  
  neg_beads <- paste0("/", beads_type, "/neg")
  neg_beads_mfi <- f.pop_MFI(gs_pop_get_data(gs[[grep("Unstained", sampleNames(gs))]], neg_beads)[[1]])
  
  
  pos_beads_mfi_substracted <- lapply(pos_beads_mfi, function(x) {
    x - neg_beads_mfi
  })
  
  pos_beads_mfi_normalized <- lapply(pos_beads_mfi_substracted, function(x) {
    b <- x[grep(" ", names(x))]
    x / b
  })
  
  pos_beads_mfi_normalized <- do.call(rbind, pos_beads_mfi_normalized)
  pos_beads_mfi_sub <- do.call(rbind, pos_beads_mfi_substracted)
  
  rownames(pos_beads_mfi_normalized) <- StrExtract(rownames(pos_beads_mfi_normalized), " ", 3)
  
  rownames(pos_beads_mfi_sub) <- StrExtract(rownames(pos_beads_mfi_sub), " ", 3)
  
  colnames(pos_beads_mfi_normalized) <- StrExtract(colnames(pos_beads_mfi_normalized), " ", 1)
  #colnames(pos_beads_mfi_normalized) <- gsub("-A", "", colnames(pos_beads_mfi_normalized))
  rownames(pos_beads_mfi_normalized) <- paste0(rownames(pos_beads_mfi_normalized), "-A")
  
  
  if(spectral == F) {
    pos_beads_mfi_normalized <- pos_beads_mfi_normalized[, colnames(pos_beads_mfi_normalized) %in% rownames(pos_beads_mfi_normalized)]
  }else {
    pos_beads_mfi_normalized <- pos_beads_mfi_normalized[, grepl("-A", colnames(pos_beads_mfi_normalized))]
  }
  list(matrix = pos_beads_mfi_normalized, 
       MFImatrix = pos_beads_mfi_sub,
       event_counts = event_conts)
}


library(lattice)
library(latticeExtra)
f.lattice_heatmap <- function(data, xlab = "", ylab = "", main = "Heatmap", colorkey = TRUE) {
  # Reverse the order of columns to reverse the x-axis
  data <- t(data) #[, ncol(data):1]
  data <- data[, ncol(data):1]
  # Define group boundaries for rows (vertical lines)
  row_group_boundaries <- c(
    UV = sum(startsWith(rownames(data), "UV")),
    V = sum(startsWith(rownames(data), "V")),
    B = sum(startsWith(rownames(data), "B")),
    YG = sum(startsWith(rownames(data), "YG")),
    R = sum(startsWith(rownames(data), "R"))
  )
  row_group_boundaries <- cumsum(row_group_boundaries) + 0.5
  
  # Define group boundaries for columns (horizontal lines)
  # Since columns are reversed, we need to adjust the positions
  col_group_boundaries <- c(
    UV = sum(startsWith(colnames(data), "UV")),
    V = sum(startsWith(colnames(data), "V")),
    B = sum(startsWith(colnames(data), "B")),
    YG = sum(startsWith(colnames(data), "YG")),
    R = sum(startsWith(colnames(data), "R"))
  )
  col_group_boundaries <- cumsum(col_group_boundaries) + 0.5
  col_group_boundaries <- ncol(data) - col_group_boundaries + 1  # Adjust for reversed columns
  
  # Define the color scale with white at 0
  # Find the maximum absolute value in the data to set symmetric limits
  max_abs_value <- max(abs(data), na.rm = TRUE)
  color_breaks <- seq(-max_abs_value, max_abs_value, length.out = 100)
  
  # Create the heatmap using levelplot
  heatmap_plot <- levelplot(
    data,
    xlab = xlab,  # Set x-axis label (empty by default)
    ylab = ylab,  # Set y-axis label (empty by default)
    main = main,  # Set the title of the heatmap
    colorkey = colorkey,  # Show or hide the color key
    scales = list(
      x = list(rot = 90, alternating = 1),  # Rotate x-axis tick labels by 90 degrees
      y = list(alternating = 1)  # Ensure consistent y-axis appearance
    ),
    par.settings = list(
     # axis.line = list(col = "transparent"),  # Remove axis lines
      axis.text = list(cex = 0.6)  # Adjust font size of axis tick labels
    ),
    at = color_breaks,  # Define the color breaks
    col.regions = colorRampPalette(c("blue", "white", "red"))(100)  # Diverging color scheme
  )
  
  # Add group boundaries
  heatmap_plot + layer(
    panel.abline(
      v = row_group_boundaries[-length(row_group_boundaries)],  # Vertical lines for row groups
      h = col_group_boundaries[-length(col_group_boundaries)],    # Horizontal lines for column groups (adjusted)
      col = "black", 
      lwd = 1
    ),
    data = list(
      row_group_boundaries = row_group_boundaries,
      col_group_boundaries = col_group_boundaries
    )  # Pass boundaries to the layer
  )
}


f.base_heatmap <- function(data, xlab = "", ylab = "", main = "Heatmap", colorkey = TRUE) {
  # Reverse the order of columns to reverse the x-axis
  data <- t(data)
  data <- data[, ncol(data):1]
  
  # Define group boundaries for rows (vertical lines)
  row_group_boundaries <- c(
    UV = sum(startsWith(rownames(data), "UV")),
    V = sum(startsWith(rownames(data), "V")),
    B = sum(startsWith(rownames(data), "B")),
    YG = sum(startsWith(rownames(data), "YG")),
    R = sum(startsWith(rownames(data), "R"))
  )
  row_group_boundaries <- cumsum(row_group_boundaries) + 0.5
  
  # Define group boundaries for columns (horizontal lines)
  # Since columns are reversed, we need to adjust the positions
  col_group_boundaries <- c(
    UV = sum(startsWith(colnames(data), "UV")),
    V = sum(startsWith(colnames(data), "V")),
    B = sum(startsWith(colnames(data), "B")),
    YG = sum(startsWith(colnames(data), "YG")),
    R = sum(startsWith(colnames(data), "R"))
  )
  col_group_boundaries <- cumsum(col_group_boundaries) + 0.5
  col_group_boundaries <- ncol(data) - col_group_boundaries + 1  # Adjust for reversed columns
  
  # Define the color scale with white at 0
  max_abs_value <- max(abs(data), na.rm = TRUE)
  color_breaks <- seq(-max_abs_value, max_abs_value, length.out = 101)  # 101 breaks for 100 colors
  color_palette <- colorRampPalette(c("blue", "white", "red"))(100)  # 100 colors
  
  # Plot the heatmap
  par(mar = c(5, 5, 4, 2))  # Adjust margins for labels
  
  image(
    1:nrow(data), 1:ncol(data), data,
    col = color_palette,
    breaks = color_breaks,  # Use breaks with one more element than colors
    xlab = xlab, ylab = ylab, main = main,
    axes = FALSE
  )
  
  # Add custom axes on all four sides
  axis(1, at = 1:nrow(data), labels = rownames(data), las = 2, cex.axis = 0.6)  # Bottom x-axis
  axis(2, at = 1:ncol(data), labels = colnames(data), las = 1, cex.axis = 0.6)  # Left y-axis
  axis(3, at = 1:nrow(data), labels = rownames(data), las = 2, cex.axis = 0.6)  # Top x-axis
  axis(4, at = 1:ncol(data), labels = colnames(data), las = 1, cex.axis = 0.6)  # Right y-axis
  
  # Add group boundaries
  abline(v = row_group_boundaries[-length(row_group_boundaries)], col = "black", lwd = 1)
  abline(h = col_group_boundaries[-length(col_group_boundaries)], col = "black", lwd = 1)
  
  # Add color key (legend)
  if (colorkey) {
    legend_image <- as.raster(matrix(color_palette, nrow = 1))
    rasterImage(legend_image, xleft = nrow(data) + 1, ybottom = 1, xright = nrow(data) + 2, ytop = ncol(data))
    text(x = nrow(data) + 2.5, y = seq(1, ncol(data), length.out = 5),
         labels = round(seq(-max_abs_value, max_abs_value, length.out = 5), 2),
         cex = 0.6, pos = 4)
  }
}







#' Extract PMT Voltage Information from a Flow Cytometry Dataset
#'
#' This function retrieves the voltage values for each photomultiplier tube (PMT) detector
#' in a flow cytometry dataset. It searches for voltage information stored in the standard
#' `$P[i]V` keyword format within the flow cytometry metadata. If a voltage keyword is not
#' found for a specific PMT, the function returns `NA` for that detector.
#'
#' @param fs A flow cytometry dataset (e.g., a `flowSet` or similar object) containing
#'           the data and metadata (keywords) for the experiment.
#'
#' @return A named numeric vector where:
#' \itemize{
#'   \item **Names**: Correspond to the names of the PMT detectors (e.g., fluorescence channels).
#'   \item **Values**: Represent the voltage values for each PMT. If a voltage value is not found
#'         for a specific PMT, the value is set to `NA`.
#' }
#'
#' @examples
#' \dontrun{
#' # Load a flow cytometry dataset
#' library(flowCore)
#' fs <- read.flowSet("path_to_fcs_files")
#'
#' # Extract PMT voltages
#' pmt_voltages <- f.PMT(fs)
#'
#' # Print the results
#' print(pmt_voltages)
#' }
#'
#' @export

f.PMT <- function(fs) {
  pdata <- parameters(fs)
  pmt_names <- pdata@data$name
  voltage_info <- sapply(1:length(pmt_names), function(i) {
    keyword <- paste0("$P", i, "V") # Construct the keyword for voltage
    if (keyword %in% names(keyword(fs))) {
      as.numeric(keyword(fs)[[keyword]])
    } else {
      return(NA) # Return NA if the keyword is not found
    }
  })
  names(voltage_info) <- pmt_names
  voltage_info
}















#' Add new channels to a flowFrame object
#'
#' @param input_fcs A flowFrame object to which the new channels will be added.
#' @param data_to_add A data frame or matrix containing the new channels to be added.
#' @param channel_names A character vector of names for the new channels.
#' @param rescale_data Logical, whether to rescale the new channels to the range [0, 100]. Default is TRUE.

#'
#' @return The updated flowFrame object with the new channels added.
f.addChannel <- function(input_fcs, data_to_add, channel_names, rescale_data = F) {
  # Ensure the number of cells matches
  stopifnot(dim(input_fcs)[1] == nrow(data_to_add))

  if (missing(channel_names)) {
    channel_names <- colnames(data_to_add)
  }

  # Rescale the data if requested
  if (rescale_data == T) {
    data_to_add <- apply(data_to_add, 2, function(x) scales::rescale(x, c(0, 1024)))
  }

  # Get the parameters of the input flowFrame
  params <- parameters(input_fcs)
  pd <- pData(params)

  for (h in 1:ncol(data_to_add)) {
    # Add a new channel
    dat <- data_to_add[, h]
    # dat <- na.omit(data_to_add[, h])
    range <- round(range(dat), 1)

    channel_number <- ncol(input_fcs) + 1
    channel_id <- paste("$P", channel_number, sep = "")
    channel_name <- channel_names[h]
    channel_range <- range
    channel_desc <- channel_name
    plist <- list(channel_name, channel_desc, range[2], range[1], range[2])

    pdif <- length(pd) - length(plist)
    # if(pdif != 0) {plist <- c(plist, rep(NA, pdif))

    # plist <- data.frame(t(plist))
    if (pdif != 0) {
      plist <- c(plist, rep(list(NA), pdif))
    }

    names(plist) <- colnames(pd)
    pd <- rbind(pd, plist)
    rownames(pd)[nrow(pd)] <- channel_id
    pData(params) <- pd

    # Update the data
    d_data <- cbind(exprs(input_fcs), data_to_add[, h, drop = F])
    input_fcs <- flowFrame(d_data, params) # , description = description(input_fcs) )

    # Add additional metadata
    keyval <- list()
    keyval[[paste("$P", channel_number, "B", sep = "")]] <- "32" # Byte
    keyval[[paste("$P", channel_number, "R", sep = "")]] <- toString(channel_range) # channel range
    keyval[[paste("$P", channel_number, "N", sep = "")]] <- channel_name
    keyval[[paste("$P", channel_number, "S", sep = "")]] <- paste0(channel_name)
    keyval[[paste("$P", channel_number, "G", sep = "")]] <- "1"
    keyval[[paste("$P", channel_number, "E", sep = "")]] <- "0,0" # amplification
    keyval[[paste("$P", channel_number, "DISPLAY", sep = "")]] <- "LIN"
    keyword(input_fcs) <- c(keyword(input_fcs), keyval)
  }

  return(input_fcs)
}

#' Add PCA-derived channels to a flowFrame object
#'
#' @param obj The input flowFrame object to which the new PCA-derived channels will be added.
#' @param pca.markers A vector of column names (or indices) of the markers in the flowFrame object that should be used to perform the PCA.
#' @param pca_loadings (Optional) The loadings (rotation matrix) from a previous PCA. If provided, the function will use these loadings to project the data onto the principal components.
#' @param return_PCs The principal components to be added as new channels to the flowFrame object. Default is 1:3 (the first three principal components).
#' @param center_data Logical, whether to center the data before performing PCA. Default is TRUE.
#' @param scale_data Logical, whether to scale the data before performing PCA. Default is FALSE, as it's usually more appropriate to compare between samples without scaling.
#' @param rescale_pca_channel Logical, whether to rescale the new PCA-derived channels to the range [0, 100]. Default is FALSE.
#'
#' @return A list with two elements:
#' - flowframe: The updated flowFrame object with the new PCA-derived channels added.
#' - loadings: The rotation matrix (loadings) from the PCA, if it was performed internally. If pca_loadings was provided, this element will be NULL.
f.addPCAchannel <- function(obj,
                            pca.markers,
                            pca_loadings,
                            return_PCs = 1:3,
                            center_data = T,
                            scale_data = F, # Default F so its faire to compare between samples
                            rescale_pca_channel = F) {
  # Check if pca_loadings is provided
  if (missing(pca_loadings)) {
    # If pca_loadings is not provided, perform PCA on the specified markers

    if (class(obj)[1] == "matrix") {
      dat <- obj[, pca.markers]

      dat_pca <- prcomp(dat, center = center_data, scale. = scale_data)
      # Project the data onto the first few principal components (specified by return_PCs)
      pca1to3_proj <- dat_pca$x[, return_PCs]
      # Add the PCA channels to the flowFrame object
      obj <- NULL
      # Store the rotation matrix (loadings)
      rotation <- dat_pca$rotation
    } else if (class(obj)[1] == "flowFrame") {
      dat <- exprs(obj)[, pca.markers]

      dat_pca <- prcomp(dat, center = center_data, scale. = scale_data)
      # Project the data onto the first few principal components (specified by return_PCs)
      pca1to3_proj <- dat_pca$x[, return_PCs]
      # Add the PCA channels to the flowFrame object
      obj <- f.addChannel(obj, pca1to3_proj, rescale_data = rescale_pca_channel)
      # Store the rotation matrix (loadings)
      rotation <- dat_pca$rotation
    }
  } else {
    # If pca_loadings is provided, use it to project the data onto the principal components
    dat <- exprs(obj)[, rownames(pca_loadings)]
    # Ensure the data and loadings have matching column names
    stopifnot(identical(colnames(dat), rownames(pca_loadings)))
    pca <- dat %*% pca_loadings
    # Add the PCA channels to the flowFrame object
    obj <- f.addChannel(obj, pca[, return_PCs], rescale_data = rescale_pca_channel)
    # Set the rotation to NULL, as it's not needed in this case
    rotation <- NULL
  }

  # Return the updated flowFrame object and the loadings (if calculated)
  return(list(flowframe = obj, loadings = rotation))
}
















#' @param fs a flow set object containing flow cytometry data
#' @param channel an optional parameter specifying the channel(s) to analyze
#' @param truncate a boolean indicating whether to truncate values above 262144 (default is true)
f.remove_tooBright_events <- function(fs, channel, truncate = T) {
  # Extract the expression data from the flow set object
  dat <- exprs(fs)

  # Find the indices of the "SC-" and "ime" channels in the data
  l_channels <- grep("SC-|ime", colnames(dat))

  # If no channel is specified, analyze all channels except the "SC-" and "ime" channels
  if (missing(channel)) {
    # Apply a function to each row of the data (excluding the "SC-" and "ime" channels)
    # The function counts the number of values greater than 262144 in each row
    outliers <- apply(dat[, -l_channels], 1, function(x) {
      sum(x > 262144)
    })

    # Find the rows where there are no values greater than 262144
    in_scale <- which(outliers == 0)

    # Subset the flow set object to only include the rows that are in scale
    fs <- fs[in_scale, ]
  }
  # If a channel is specified, analyze that channel
  else {
    # Apply the same function as above, but only to the specified channel
    outliers <- apply(dat[, channel, drop = F], 1, function(x) {
      sum(x > 262144)
    })

    # Find the rows where there are no values greater than 262144
    in_scale <- which(outliers == 0)

    # Subset the flow set object to only include the rows that are in scale
    fs <- fs[in_scale, ]

    # If truncate is set to true, truncate any values greater than 262144 in the non-channel columns
    if (truncate == T) {
      exprs(fs)[, -l_channels][exprs(fs)[, -l_channels] > 262144] <- 262144
    }
  }

  # Return the modified flow set object
  return(fs)
}



#' Finds markers in the FCS file

#' @param  frame: a flowframe
#' @param marker.list: a vector of markers

Find.markers <- function(frame, marker.list) {
  # Parameters:
  #* frame: a flowFrame in the flowSet
  #** marker.list: A vector of characters
  # Output:
  #* channels.ind: a vector of channels numbers in the frame  corresponding to marker.list
  if (class(frame) == "cytoframe") {
    frame <- cytoframe_to_flowFrame(frame)
  } # flowWorkspace 4 uses cytoframe instead of flwoframe, so need to convert

  channels.ind <- unlist(lapply(marker.list, function(x) {
    ind <- grep(x, frame@parameters@data[, 2], ignore.case = T)
    ind_store <- ind
    if (length(ind) == 0) {
      warning(paste(x, "not found, check markers!"))
      return(NA)
    } else {
      if (length(ind) > 1) {
        cnt <- 0
        repeat{
          cnt <- cnt + 1
          fs.markers <- unlist(lapply(frame@parameters@data[, 2], function(x) unlist(strsplit(x, " "))[cnt]))
          ind <- match(x, fs.markers)
          if (is.na(ind)) {
            fs.markers <- unlist(lapply(frame@parameters@data[, 2], function(x) unlist(strsplit(x, "-"))[cnt]))
            ind <- match(x, fs.markers)
            if (!is.na(ind)) {
              break
            }
          } else {
            break
          }
          if (cnt >= 10) {
            if (length(ind_store) >= 2) {
              ind <- grep(paste0(x, "_"), frame@parameters@data[, 2], ignore.case = T)
              if (length(ind) == 0) {
                ind <- ind_store[1]
                warning(paste(x, "found more than one, choosing first. Check markers!"))
              }
            } else {
              warning(paste(x, "not found, check markers!"))
            }
            break
          }
        }
      }
    }
    return(ind)
  }))
  names(channels.ind) <- marker.list
  # Removing NAs in the channel vector, as not all centres have all of these markers
  # Note that most centres should have Live/CD4/CD8/CD44/CD62/CD5/CD161/CD25 in their channels
  ind <- which(is.na(channels.ind))
  if (length(ind) != 0) {
    channels.ind <- channels.ind[-ind]
  }
  return(channels.ind)
}


#' Perform Logicle transformation with known w.
#' a csv file' path was supplid with a set of w for the transfromation.
#'
#' @param obj flowframe or GatingSet
#' @param trans.chans NULL by defual, transform all log channels
#' @param PathtoParameters Default NULL, a defult value 0.8 will be used.
#' Otherwise, w specified in the csv the PathtoParameters point to will be used
#' @return same object as the input with transformed data
#' PathtoParameters <- paste0(AUTOFLOW_destination, "Reference_Data/A5_trans_parameters_T.csv")
#' updated Oct 31, 2023
#' obj = fs.raw[[1]]; trans.chans = NULL; PathtoParameters = NULL; w_default = 0.8
transform.logicle <- function(obj,
                              trans.chans = NULL,
                              PathtoParameters = NULL,
                              w_default = 0.8) {
  # handling obj
  if (class(obj) == "flowFrame") {
    f <- obj
    # f.set <- f
  } else { # handling sets
    # temp <- obj
    if (class(obj) == "GatingSet") {
      # f.set <- gs_pop_get_data( obj, tail(getNodes(obj), 1) ) # get the last node
      f.set <- gs_pop_get_data(obj, "root") # get the root
    }
    if (class(obj) == "flowSet") {
      f.set <- obj
    }

    f <- f.set[[1]] # a flowset
  }

  # This line is added to fix the column selection problem introduced in setting DISPLAY
  # parameter in "LIN" on DIVA by accident
  # f@description$P35DISPLAY <- "LOG"

  if (is.null(trans.chans)) {
    # This section is changed to fix the column selection problem introduced in keep.col

    log.channels <- names(f@description)[grepl("DISPLAY", names(f@description))]

    temp_order <- gsub("([0-9]+).*$", "\\1", log.channels)
    temp_order <- as.numeric(gsub("P", "", temp_order))
    log.channels <- log.channels[order(temp_order)]

    trans.chans <- which(f@description[log.channels] == "LOG")
    trans.chans <- gsub("DISPLAY", "N", names(trans.chans))
    trans.chans <- paste("$", trans.chans, sep = "")
    trans.chans <- match(unlist(f@description[trans.chans]), colnames(f)) # Get Index
    # make sure the index don't mess up with
    if (length(trans.chans) == 0) {
      warnings("Couldn't find Log channels,
                all channels except FSC, SSC,
                and time will be transformed.")
      trans.chans <- 1:ncol(f)
      trans.chans <- trans.chans[-c(
        grep(colnames(f), pattern = "FSC*"),
        grep(colnames(f), pattern = "SSC*"),
        grep(colnames(f), pattern = "Time*")
      )]
    }
  }

  # Construct transformationList


  if (is.null(PathtoParameters)) {
    trans.list <- lapply(colnames(f)[trans.chans], function(x) {
      if (class(obj) == "flowFrame") {
        logicleTransform(w = w_default, t = 262143, m = 4.5)
      } else {
        logicle_trans(w = w_default, t = 262143, m = 4.5, a = 0)
      }
    })

    names(trans.list) <- colnames(f)[trans.chans]
  } else {
    para <- read.csv(PathtoParameters, row.names = 1)

    trans.list <- lapply(rownames(para), function(x) {
      if (class(obj) == "flowFrame") {
        logicleTransform(w = para[x, "w"], t = para[x, "t"], m = para[x, "m"], a = 0)
      } else {
        logicle_trans(w = para[x, "w"], t = para[x, "t"], m = para[x, "m"], a = 0)
      }
    })

    names(trans.list) <- rownames(para)
  }


  if (!identical(names(trans.list), colnames(f)[trans.chans])) {
    para <- para[order(match(names(trans.list), colnames(f)[trans.chans])), ]
    trans.list <- trans.list[order(match(names(trans.list), colnames(f)[trans.chans]))]
  }

  stopifnot(identical(names(trans.list), colnames(f)[trans.chans]))



  if (class(obj) == "flowFrame") {
    dat <- exprs(obj)
    n_of_na <- apply(dat, 1, function(x) {
      sum(is.na(x))
    })
    dat <- dat[n_of_na == 0, ]
    exprs(obj) <- dat




    trans <- transformList(rownames(para), trans.list) # transformerlist
    obj <- transform(obj, trans)
    obj@parameters@data$trans <- para$w[match(obj@parameters@data$name, rownames(para))]
    obj@parameters@varMetadata <- rbind(obj@parameters@varMetadata, trans = "Transformation Parameter")
  } else {
    trans <- transformerList(rownames(para), trans.list) # transformerlist
    obj <- transform(obj, trans)

    fs <- fsApply(gs_pop_get_data(obj), function(x) {
      x@parameters@data$trans <- para$w[match(x@parameters@data$name, rownames(para))]
      x@parameters@varMetadata <- rbind(x@parameters@varMetadata, trans = "Transformation Parameter")
      x
    })

    obj <- GatingSet(fs)
    # obj@transformation <- trans
  }


  return(obj)
}








#' Perform Logicle transformation with known w.
#' a csv file' path was supplid with a set of w for the transfromation.
#'
#' @param obj flowframe or GatingSet
#' @param trans.chans NULL by defual, transform all log channels
#' @param PathtoParameters Default NULL, a defult value 0.8 will be used.
#' Otherwise, w specified in the csv the PathtoParameters point to will be used
#' @return same object as the input with transformed data
#' PathtoParameters <- paste0(AUTOFLOW_destination, "Reference_Data/A5_trans_parameters_T.csv")
#' updated Oct 31, 2023
#' obj = fs.raw[[1]]; trans.chans = NULL; PathtoParameters = NULL; w_default = 0.8
#' 
#' 
#' this version is for R4, which prevent you add a column to the parameter 
transform.logicle_R4 <- function(obj,
                              trans.chans = NULL,
                              PathtoParameters = NULL,
                              w_default = 0.8) {
  # handling obj
  if (class(obj) == "flowFrame") {
    f <- obj
    # f.set <- f
  } else { # handling sets
    # temp <- obj
    if (class(obj) == "GatingSet") {
      # f.set <- gs_pop_get_data( obj, tail(getNodes(obj), 1) ) # get the last node
      f.set <- gs_pop_get_data(obj, "root") # get the root
    }
    if (class(obj) == "flowSet") {
      f.set <- obj
    }
    
    f <- f.set[[1]] # a flowset
  }
  
  # This line is added to fix the column selection problem introduced in setting DISPLAY
  # parameter in "LIN" on DIVA by accident
  # f@description$P35DISPLAY <- "LOG"
  
  if (is.null(trans.chans)) {
    # This section is changed to fix the column selection problem introduced in keep.col
    
    log.channels <- names(f@description)[grepl("DISPLAY", names(f@description))]
    
    temp_order <- gsub("([0-9]+).*$", "\\1", log.channels)
    temp_order <- as.numeric(gsub("P", "", temp_order))
    log.channels <- log.channels[order(temp_order)]
    
    trans.chans <- which(f@description[log.channels] == "LOG")
    trans.chans <- gsub("DISPLAY", "N", names(trans.chans))
    trans.chans <- paste("$", trans.chans, sep = "")
    trans.chans <- match(unlist(f@description[trans.chans]), colnames(f)) # Get Index
    # make sure the index don't mess up with
    if (length(trans.chans) == 0) {
      warnings("Couldn't find Log channels,
                all channels except FSC, SSC,
                and time will be transformed.")
      trans.chans <- 1:ncol(f)
      trans.chans <- trans.chans[-c(
        grep(colnames(f), pattern = "FSC*"),
        grep(colnames(f), pattern = "SSC*"),
        grep(colnames(f), pattern = "Time*")
      )]
    }
  }
  
  # Construct transformationList
  
  
  if (is.null(PathtoParameters)) {
    trans.list <- lapply(colnames(f)[trans.chans], function(x) {
      if (class(obj) == "flowFrame") {
        logicleTransform(w = w_default, t = 262143, m = 4.5)
      } else {
        logicle_trans(w = w_default, t = 262143, m = 4.5, a = 0)
      }
    })
    
    names(trans.list) <- colnames(f)[trans.chans]
  } else {
    para <- read.csv(PathtoParameters, row.names = 1)
    
    trans.list <- lapply(rownames(para), function(x) {
      if (class(obj) == "flowFrame") {
        logicleTransform(w = para[x, "w"], t = para[x, "t"], m = para[x, "m"], a = 0)
      } else {
        logicle_trans(w = para[x, "w"], t = para[x, "t"], m = para[x, "m"], a = 0)
      }
    })
    
    names(trans.list) <- rownames(para)
  }
  
  
  if (!identical(names(trans.list), colnames(f)[trans.chans])) {
    trans.list <- trans.list[order(match(names(trans.list), colnames(f)[trans.chans]))]
  }
  
  stopifnot(identical(names(trans.list), colnames(f)[trans.chans]))
  
  
  
  if (class(obj) == "flowFrame") {
    dat <- exprs(obj)
    n_of_na <- apply(dat, 1, function(x) { 
      sum(is.na(x))
    })
    dat <- dat[n_of_na == 0, ]
    exprs(obj) <- dat
    
    
    
    
    trans <- transformList(names(trans.list), trans.list) # transformerlist
    obj <- transform(obj, trans)

  } else {
    trans <- transformerList(names(trans.list), trans.list) # transformerlist
    obj <- transform(obj, trans)
    
    # obj@transformation <- trans
  }
  
  
  return(obj)
}



# Margin removal function

removeMargins <- function(f, chans, sens = 1, debris = FALSE, neg = 500, verbose = T, return.gate = F) {
  neg <- rep(neg, length(chans))
  data <- exprs(f)
  margins <- c()
  marg.gates <- c()
  if (!debris) {
    for (chan in chans)
    {
      if (is.character(chan)) {
        chan <- which(colnames(f) == chan)
      }
      stain.max <- max(data[, chan], na.rm = T)
      margins <- which(data[, chan] >= stain.max * sens)
      data <- data[-margins, ]
      if (verbose == T) {
        print(paste(length(margins), "margin events in", colnames(f)[chan], "will be removed.", sep = " "))
      }
      marg.gates <- append(marg.gates, stain.max * sens - 1)
    }
    if (return.gate) {
      return(marg.gates)
    }
  } else {
    for (chan in chans)
    {
      if (is.character(chan)) {
        chan <- which(colnames(f) == chan)
      }
      stain.min <- min(data[, chan], na.rm = T)
      margins <- which(data[, chan] <= stain.min * sens)
      data <- data[-margins, ]
      if (verbose == T) {
        print(paste(length(margins), "debris events in", colnames(f)[chan], "will be removed.", sep = " "))
      }
    }
  }
  for (i in 1:length(chans))
  {
    if (neg[i] < 500) {
      ch <- ifelse(is.character(chans[i]), yes = which(colnames(f) == chans[i]), no = chans[i])
      negs <- which(data[, ch] < neg[i])
      margins <- negs
      if (length(margins) != 0) {
        data <- data[-margins, ]
      }
      if (verbose == T) {
        print(paste(length(margins), "negative events in", colnames(f)[ch], "will be removed.", sep = " "))
      }
    }
  }
  exprs(f) <- data

  return(f)
}


# Data rotation

rotate.data <- function(data, chans = NULL, theta = NULL, min.max = F) {
  if (nrow(data) < 3) {
    print("Cannot rotate a matrix with less than 3 rows, returning back the original matrix")
    return(list(data = data, theta = theta))
  } else {
    if (inherits(data, "flowFrame") & !is.null(chans)) {
      if (all(is.na(exprs(data)[, 1]))) {
        return("Cannot rotate a flowFrame with all cells NA")
      }
      no.nas <- which(!is.na(exprs(data)[, chans[1]]))
      data.new <- exprs(data)[no.nas, chans]
      if (is.null(theta)) {
        reg.slope <- atan(lm(data.new[, 1] ~ data.new[, 2])$coefficients[2])
        slope <- atan((max(data.new[, 2]) - min(data.new[, 2])) / (max(data.new[, 1]) - min(data.new[, 1])))
        theta <- ifelse(min.max, no = pi / 2 - reg.slope, yes = slope - pi / 2)
      }
      data.new <- data.new %*% matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2, byrow = T)
      exprs(data)[no.nas, chans] <- data.new
    } else {
      col.names <- colnames(data)
      data <- data %*% matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2, byrow = T)
      colnames(data) <- col.names
    }
    return(list(data = data, theta = theta))
  }
}


# flowDensity rotation


# Rotating back a flowDensity object
rotate.fd <- function(fd.object, angle) {
  new.f <- new(Class = "CellPopulation")
  no.na <- which(!is.na(exprs(fd.object@flow.frame)[, 1]))
  dat <- fd.object@flow.frame
  temp <- rotate.data(dat[no.na, ], fd.object@channels, theta = -angle)$data
  exprs(dat)[no.na, ] <- exprs(temp)
  new.f@flow.frame <- dat
  new.f@filter <- rotate.data(fd.object@filter, theta = -angle)$data
  colnames(new.f@filter) <- colnames(fd.object@filter)
  new.f@channels <- fd.object@channels
  new.f@proportion <- fd.object@proportion
  new.f@cell.count <- fd.object@cell.count
  new.f@index <- fd.object@index
  return(new.f)
}


f.compare_comp_matrix <- function(fs, comp) {
  machine_comp <- spillover(fs)
  machine_comp <- machine_comp[!unlist(lapply(machine_comp, is.null))][[1]]



  comp <- round(comp * 100, 2)
  machine_comp <- round(machine_comp * 100, 2)
  comp_diff <- comp - machine_comp
  # comp_diff <- comp_diff[, apply(comp_diff,2, sum)!=0 ]

  notes <- lapply(seq(ncol(comp_diff)), function(i) {
    receiver <- colnames(comp_diff)[i]
    emitters_index <- which(comp_diff[, i] != 0)
    emitters <- rownames(comp_diff)[emitters_index]

    do.call(c, lapply(seq(emitters), function(j) {
      paste0(
        "", emitters[j], " to ", receiver, ": ",
        machine_comp[emitters_index[j], i], " to ", comp[emitters_index[j], i], "  "
      )
    }))
  })
  notes <- do.call(c, notes)

  if (is.null(notes)) {
    notes <- "No COMP tweak applied"
  }
  notes
}




# Function that draw flat side gate
# x a two column boundary
f.NiceGate <- function(x, top = T, left = T, topline = 4.0, leftline = 0, rightline =4, bottomline = 0) {
  min.x.index <- which.min(x[, 1])
  min.y.index <- which.min(x[, 2])
  max.x.index <- which.max(x[, 1])
  max.y.index <- which.max(x[, 2])

  if (is.na(topline)) {
    top_value <- x[max.y.index, 2]
  } else {
    top_value <- max(topline, x[max.y.index, 2])
  }

  if (is.na(leftline)) {
    leftvalue <- x[min.x.index, 1]
  } else {
    leftvalue <- min(leftline, x[min.x.index, 1])
  }

  if (is.na(rightline)) {
    rightvalue <- x[min.x.index, 1]
  } else {
    rightvalue <- min(rightline, x[min.x.index, 1])
  }

  if (is.na(bottomline)) {
    bottomvalue <- x[min.x.index, 1]
  } else {
    bottomvalue <- min(bottomline, x[min.x.index, 1])
  }
  
  

  if (top == T) {
    left.up.corner <- c(leftvalue, top_value)
    right.up.corner <- c(x[1, 1], top_value)

    if (left == T) {
      left.down.corner <- c(leftvalue, x[min.y.index, 2])
      # right.corner <- c(x[1, 1], x[max.y.index, 2])
      x <- rbind(x[1:min.y.index, ], left.down.corner, left.up.corner, right.up.corner, x[1, ])
    } else {
      x <- rbind(x[1:min.x.index, ], left.up.corner, right.up.corner, x[1, ])
    }
  }


  x
}





#' compensation by defaul or using a dedicated csv file



f.custom_compensation <- function(fs, csv) {
  if (!is.null(csv)) {
    comp <- tryCatch(
      {
        comp <- csv
      },
      error = function(e) {
        return(NULL)
      }
    )
  } else {
    comp <- NULL
  }

  machine_comp <- spillover(fs)
  machine_comp <- machine_comp[!unlist(lapply(machine_comp, is.null))][[1]]

  if (is.null(comp)) {
    comp <- machine_comp
    rownames(comp) <- colnames(comp)

    comp_by_pct <- F
    if (comp_by_pct == T) {

    }
  } else {
    if (is.null(rownames(comp)[[1]])) {
      rownames(comp) <- colnames(comp)
    } else {
      rownames(comp) <- StrExtract(rownames(comp), " :: ", 1)
      colnames(comp) <- rownames(comp) # replace the colnames since colnames were messed up read from the flowjo direct exportef file
    }


    if (identical(sort(colnames(comp)), sort(colnames(machine_comp)))) { # Make sure the fluorochomes are the same
      if (identical(rownames(comp), colnames(machine_comp)) == F) {
        comp <- comp[colnames(machine_comp), colnames(machine_comp)] # Re order if the order is diferent
      }
    }
  }

  return(comp)
}

#' x <- fs[[1]]; show.plots = T

f.singlets_gating <- function(x, show.plots = T) {
  # flowPlot(x, channels = c("FSC-A", "FSC-H"))
  centered_data <- scale(exprs(x)[, c("FSC-A", "FSC-H")], scale = FALSE) # Center the data (subtract means)
  slope <- coef(lm(centered_data[, 2] ~ centered_data[, 1]))[2] # Slope of the data
  rot.theta <- 0.7 # atan(slope) *1.3  # Angle to rotate by

  # rotate the data
  tp_fs.rot_HA <- rotate.data(x, c("FSC-A", "FSC-H"), theta = rot.theta)$data
  # flowPlot(tp_fs.rot_HA, channels = c("FSC-A", "FSC-H"))
  # plotDens(tp_fs.rot_HA, c("FSC-A", "FSC-H"),density.overlay = c(T,T))

  # find the cutting point for x axis
  #pk <- getPeaks(tp_fs.rot_HA, "FSC-A", tinypeak.removal = 1 / 10)

    tp_ct_FSCA <- deGate(tp_fs.rot_HA, "FSC-A", tinypeak.removal = 1 / 10, use.upper = T, upper = F)

    tp_ct_FSCA <- c(tp_ct_FSCA, deGate(tp_fs.rot_HA, "FSC-A", tinypeak.removal = 1 / 10, all.cuts = T))
    

  tp_ct_FSCA <- min(tp_ct_FSCA)


  # plotDens(tp_fs.rot_HA, c("FSC-A", "FSC-H"),density.overlay = c(T,T))
  # find the cutting point for y axis
  # tp_pk_FSCH <- getPeaks(tp_fs.rot_HA, "FSC-H", tinypeak.removal = 1 / 10)$Peaks

  # find the up limit of FSC-H
  tp_ct_FSCH_up <- deGate(tp_fs.rot_HA, "FSC-H", tinypeak.removal = 1 / 10, use.upper = T, upper = T)
  # find the lower limit of FSC-H
  tp_ct_FSCH_lo <- deGate(tp_fs.rot_HA, "FSC-H", tinypeak.removal = 1 / 10, use.upper = T, upper = F)

  # Perform gating
  tp_cp_FSCA_FSCH_TT <- flowDensity(tp_fs.rot_HA, c("FSC-A", "FSC-H"), position = c(T, T), gates = c(tp_ct_FSCA, tp_ct_FSCH_lo))
  tp_cp_FSCA_FSCH_NAF <- flowDensity(tp_cp_FSCA_FSCH_TT@flow.frame, c("FSC-A", "FSC-H"), position = c(NA, F), gates = c(NA, tp_ct_FSCH_up))

  tp_cp_singlets <- rotate.fd(tp_cp_FSCA_FSCH_NAF, angle = rot.theta)

  if (show.plots == T) {
    par(mfrow = c(1, 4))
    flowPlot(tp_fs.rot_HA, channels = c("FSC-A", "FSC-H"), ylim = c(-50000, 20000), xlim = c(0, 300000))
    flowPlot(x, channels = c("FSC-A", "FSC-H"), xlim = c(0, 250000), ylim = c(0, 250000))
    flowPlot(tp_cp_FSCA_FSCH_NAF, channels = c("FSC-A", "FSC-H"), ylim = c(-50000, 20000), xlim = c(0, 300000))

    flowPlot(x, channels = c("FSC-A", "FSC-H"), xlim = c(0, 250000), ylim = c(0, 250000))
    lines(tp_cp_singlets@filter, col = "green", lty = 1, lwd = 2)
  }



  return(polygonGate(filterId = "Singlets", .gate = tp_cp_singlets@filter))
}

#
#' fe
#'
#'
#' debug
#' x = gs_pop_get_data(gs, "Singlets")[[4]]

f.livedead_gating <- function(x, show.plots = T) {
  tp_ct_dead <- 4 # deGate(x, "UV446-A", tinypeak.removal = 1 / 25, upper = T)
  ## flowPlot(x, channels = c("UV379-A", "UV446-A"))
  ## plotDens(x, c("UV379-A", "UV446-A"), density.overlay = c(T,T)  )

  tp_cp_live <- flowDensity(x, c("UV379-A", "UV446-A"), position = c(NA, F), gates = c(NA, tp_ct_dead))


  if (show.plots == T & nrow(x) >10) {
    par(mfrow = c(1, 4))
    markers <- c("UV379-A", "UV446-A")
    flowPlot(x, channels = markers)
    flowPlot(x, channels = markers)
    lines(tp_cp_live@filter, col = "green", lty = 1, lwd = 2)
  }

  return(polygonGate(filterId = "Live", .gate = tp_cp_live@filter))
}





#' @param x a flowrame
#' @param Tissue tissue type, if missing, the function will extract from the name of the file

# x = fs[[11]]; Tissue = "PB"
f.Expanded_Lymph_gating <- function(x, # obj,
                                    Tissue,
                                    show.plots = T) {
  # get tissue infor based on the flowframe discription data
  if (missing(Tissue)) {
    Tissue <- c("BM", "PB", "FL", "LN", "OT")

    Tissue <- Tissue[unlist(lapply(Tissue, function(t) {
      grepl(t, identifier(x))
    }))]
  }

  if (Tissue %in% c("BM", "PB", "FL")) {
    # par(mfrow = c(1, 1))
    # x <-  fs.clean[[1]]
    print(identifier(x))

    # 1. FIND BEST NLRBC GATING ON CD45
    ## 1A. log transform the SSC-A channel for better population identification
    exprs(x)[, "SSC-A"] <- log10(exprs(x)[, "SSC-A"]) # log ssc-A
    # flowPlot(x, channels = c("UV379-A", "SSC-A"))

    # 2. find a cutoff to gate wbc
    # tried a few ways, the all cut seems works best
    ct_CD45_UV379_NLRBC_preSet <- 1.2
    pks <- getPeaks(x, "UV379-A", tinypeak.removal = 1 / 25)
    if (sum(pks$Peaks < 1) == 0) {
      tp_ct_45_NLRBC <- deGate(x, "UV379-A", use.upper = T, upper = F, tinypeak.removal = 1 / 100)
    } else {
      tp_ct_45_NLRBC <- deGate(x, "UV379-A", all.cuts = T, tinypeak.removal = 1 / 100)
    }

    tp_ct_45_NLRBC <- tp_ct_45_NLRBC[which.min(abs(tp_ct_45_NLRBC - ct_CD45_UV379_NLRBC_preSet))]

    # Perform gating
    tp_cp_45_wbc <- flowDensity(x, c("UV379-A", "SSC-A"), position = c(T, NA), gates = c(tp_ct_45_NLRBC, NA))
    ## flowPlot(tp_cp_45_wbc, channels = c("UV379-A", "SSC-A"))





    # rotate gate to gate lymph and grans
    rotate.theta <- 0.68
    tp_fs.rot_45SSCA <- rotate.data(getflowFrame(tp_cp_45_wbc), c("UV379-A", "SSC-A"), theta = rotate.theta)$data
    # flowPlot(tp_fs.rot_45SSCA, c("UV379-A", "SSC-A"), density.overlay = c(T,T)  )

    # Expaned Lymph Gate
    tp_ct_SSC_rot <- deGate(tp_fs.rot_45SSCA, "SSC-A", all.cuts = T)

    ct_SSC_rot_preSet <- 2
    tp_ct_SSC_rot <- tp_ct_SSC_rot[which.min(abs(tp_ct_SSC_rot - ct_SSC_rot_preSet))]


    if (tp_ct_SSC_rot > ct_SSC_rot_preSet) { # case with a lot of grans, masking lymphs,
      # typicl cases: BF25-0177
      tp_ct_SSC_rot <- deGate(tp_fs.rot_45SSCA, "SSC-A", use.upper = F, upper = F)
    }

    tp_cp_SSC_expLymph.rot <- flowDensity(tp_fs.rot_45SSCA, c("UV379-A", "SSC-A"), position = c(NA, F), gates = c(NA, tp_ct_SSC_rot))
    # plotDens(tp_cp_SSC_expLymph.rot, c("UV379-A", "SSC-A"), density.overlay = c(T,T)  )

    tp_cp_SSC_expLymph <- rotate.fd(tp_cp_SSC_expLymph.rot, angle = rotate.theta)
    # plotDens(tp_cp_SSC_expLymph, c("UV379-A", "SSC-A"), density.overlay = c(T,T)  )

    # Grans Gate
    tp_cp_SSC_grans.rot <- flowDensity(tp_fs.rot_45SSCA, c("UV379-A", "SSC-A"), position = c(NA, T), gates = c(NA, tp_ct_SSC_rot))
    # plotDens(tp_cp_SSC_expLymph.rot, c("UV379-A", "SSC-A"), density.overlay = c(T,T)  )

    tp_cp_SSC_grans <- rotate.fd(tp_cp_SSC_grans.rot, angle = rotate.theta)

    tp_ct_SSC_grans <- deGate(tp_cp_SSC_grans, "SSC-A", use.upper = T, upper = F)
    tp_cp_SSC_grans <- flowDensity(tp_cp_SSC_grans, c("UV379-A", "SSC-A"), position = c(NA, T), gates = c(NA, tp_ct_SSC_grans))

    # plotDens(tp_cp_SSC_grans, c("UV379-A", "SSC-A"), density.overlay = c(T,T)  )


    if (show.plots == T) {
      par(mfrow = c(1, 4))
      flowPlot(x, channels = c("UV379-A", "SSC-A"), xlim = c(-1, 4.5), title = identifier(x))
      flowPlot(tp_cp_45_wbc, channels = c("UV379-A", "SSC-A"), xlim = c(-1, 4.5), title = "wbc")
      flowPlot(tp_fs.rot_45SSCA, channels = c("UV379-A", "SSC-A"))
      abline(h = tp_ct_SSC_rot)
      flowPlot(x, channels = c("UV379-A", "SSC-A"), xlim = c(-1, 4.5))
      lines(tp_cp_SSC_expLymph@filter, col = "blue", lty = 1, lwd = 2)
      lines(tp_cp_SSC_grans@filter, col = "green", lty = 1, lwd = 2)
    }

    exprs(x)[, "SSC-A"] <- 10^(exprs(x)[, "SSC-A"]) # log ssc-A
    tp_cp_SSC_expLymph@filter[, "SSC-A"] <- 10^tp_cp_SSC_expLymph@filter[, "SSC-A"]
    tp_cp_SSC_grans@filter[, "SSC-A"] <- 10^tp_cp_SSC_grans@filter[, "SSC-A"]
    # Translate the filter to gate
  } else { # LN gating

    exprs(x)[, "SSC-A"] <- log10(exprs(x)[, "SSC-A"]) # log ssc-A

    ct_CD45_UV379_NLRBC_preSet <- 1.2

    pks <- getPeaks(x, "UV379-A", tinypeak.removal = 1 / 25)
    if (sum(pks$Peaks < 1) == 0) {
      tp_ct_45_NLRBC <- deGate(x, "UV379-A", use.upper = T, upper = F, tinypeak.removal = 1 / 100)
    } else {
      tp_ct_45_NLRBC <- deGate(x, "UV379-A", all.cuts = T, tinypeak.removal = 1 / 100)
    }


    tp_ct_45_NLRBC <- tp_ct_45_NLRBC[which.min(abs(tp_ct_45_NLRBC - ct_CD45_UV379_NLRBC_preSet))]

    # Perform gating
    tp_cp_45_wbc <- flowDensity(x, c("UV379-A", "SSC-A"), position = c(T, NA), gates = c(tp_ct_45_NLRBC, NA))
    ## flowPlot(tp_cp_45_wbc, channels = c("UV379-A", "SSC-A"))


    rotate.theta <- 0.68
    ct_SSC_rot_preSet <- 2
    tp_fs.rot_45SSCA <- rotate.data(getflowFrame(tp_cp_45_wbc), c("UV379-A", "SSC-A"), theta = rotate.theta)$data
    # plotDens(tp_fs.rot_45SSCA, c("UV379-A", "SSC-A"), density.overlay = c(T,T)  )



    pks <- getPeaks(tp_fs.rot_45SSCA, "SSC-A")

    if (length(pks$Peaks) == 1) {
      if (pks$Peaks > ct_SSC_rot_preSet) {
        tp_ct_SSC_rot <- deGate(tp_fs.rot_45SSCA, "SSC-A", use.upper = T, upper = T)
      } else {
        tp_ct_SSC_rot <- deGate(tp_fs.rot_45SSCA, "SSC-A", use.upper = F, upper = T)
      }
    } else {
      tp_ct_SSC_rot <- deGate(tp_fs.rot_45SSCA, "SSC-A")
    }

    # Expaned Lymph Gate

    if (tp_ct_SSC_rot > 2.5) { # case with a lot of grans, masking lymphs,
      # typicl cases: BF25-0177
      tp_ct_SSC_rot <- deGate(tp_fs.rot_45SSCA, "SSC-A", use.upper = F, upper = F)
    }

    tp_cp_SSC_expLymph.rot <- flowDensity(tp_fs.rot_45SSCA, c("UV379-A", "SSC-A"), position = c(NA, F), gates = c(NA, tp_ct_SSC_rot))
    # plotDens(tp_cp_SSC_expLymph.rot, c("UV379-A", "SSC-A"), density.overlay = c(T,T)  )

    tp_cp_SSC_expLymph <- rotate.fd(tp_cp_SSC_expLymph.rot, angle = rotate.theta)
    # plotDens(tp_cp_SSC_expLymph, c("UV379-A", "SSC-A"), density.overlay = c(T,T)  )

    # Grans Gate
    tp_cp_SSC_grans.rot <- flowDensity(tp_fs.rot_45SSCA, c("UV379-A", "SSC-A"), position = c(NA, T), gates = c(NA, tp_ct_SSC_rot))
    # plotDens(tp_cp_SSC_expLymph.rot, c("UV379-A", "SSC-A"), density.overlay = c(T,T)  )

    tp_cp_SSC_grans <- rotate.fd(tp_cp_SSC_grans.rot, angle = rotate.theta)

    tp_ct_SSC_grans <- deGate(tp_cp_SSC_grans, "SSC-A", use.upper = T, upper = F)
    tp_cp_SSC_grans <- flowDensity(tp_cp_SSC_grans, c("UV379-A", "SSC-A"), position = c(NA, T), gates = c(NA, tp_ct_SSC_grans))

    # plotDens(tp_cp_SSC_grans, c("UV379-A", "SSC-A"), density.overlay = c(T,T)  )


    if (show.plots == T) {
      par(mfrow = c(1, 4))
      flowPlot(x, channels = c("UV379-A", "SSC-A"), xlim = c(-1, 4.5), title = identifier(x))
      flowPlot(tp_cp_45_wbc, channels = c("UV379-A", "SSC-A"), xlim = c(-1, 4.5), title = "wbc")
      flowPlot(tp_fs.rot_45SSCA, channels = c("UV379-A", "SSC-A"))
      abline(h = tp_ct_SSC_rot)
      flowPlot(x, channels = c("UV379-A", "SSC-A"), xlim = c(-1, 4.5))
      lines(tp_cp_SSC_expLymph@filter, col = "blue", lty = 1, lwd = 2)
      lines(tp_cp_SSC_grans@filter, col = "green", lty = 1, lwd = 2)
    }

    exprs(x)[, "SSC-A"] <- 10^(exprs(x)[, "SSC-A"]) # log ssc-A
    tp_cp_SSC_expLymph@filter[, "SSC-A"] <- 10^tp_cp_SSC_expLymph@filter[, "SSC-A"]
    tp_cp_SSC_grans@filter[, "SSC-A"] <- 10^tp_cp_SSC_grans@filter[, "SSC-A"]
  } # LN gating


  return(list = c(
    Exp_lymph = polygonGate(filterId = "Exp_Lymph", .gate = tp_cp_SSC_expLymph@filter),
    Grans = polygonGate(filterId = "Grans", .gate = tp_cp_SSC_grans@filter)
  ))
}



#' @param x a flowrame
f.WBC_gating <- function(x, show.plots = T) {
  
  # 1. FIND BEST NLRBC GATING ON CD45
  ## 1A. log transform the SSC-A channel for better population identification
  exprs(x)[, "SSC-A"] <- log10(exprs(x)[, "SSC-A"]) # log ssc-A
  # flowPlot(x, channels = c("UV379-A", "SSC-A"))
  
  # 2. find a cutoff to gate wbc
  # tried a few ways, the all cut seems works best
  ct_CD45_UV379_NLRBC_preSet <- 1.5
  # pks <- getPeaks(x, "UV379-A", tinypeak.removal = 1 / 25)
  # 
  # if (sum(pks$Peaks < ct_CD45_UV379_NLRBC_preSet) == 0) {
  #   tp_ct_45_NLRBC <- deGate(x, "UV379-A", use.upper = T, upper = F, tinypeak.removal = 1 / 25, alpha = 1 / 100)
  # } else {
  #   tp_ct_45_NLRBC <- deGate(x, "UV379-A", all.cuts = T, tinypeak.removal = 1 / 25)
  # }
  
  tp_ct_45_NLRBC <- f.roughGate(x, "UV379-A", thrd = 1.5, tinypeak.removal = 1/25)
  #tp_ct_45_NLRBC <- tp_ct_45_NLRBC[which.min(abs(tp_ct_45_NLRBC - ct_CD45_UV379_NLRBC_preSet))]
  
  tp_cp_45_wbc_rough <- flowDensity(x, c("UV379-A", "SSC-A"), position = c(F, NA), gates = c(tp_ct_45_NLRBC, NA))
  
  
  if (show.plots == T) {
    par(mfrow = c(1, 4))
    flowPlot(x, channels = c("UV379-A", "SSC-A"), 
             xlim = c(-1, 4.5), ylim = c(2.5, 5.5), title = identifier(x))
    abline(v = tp_ct_45_NLRBC)
    
    flowPlot(tp_cp_45_wbc_rough, channels = c("UV379-A", "SSC-A"),
             xlim = c(-1, 4.5), ylim = c(2.5, 5.5), title = identifier(x))
  }
  
  
  tp_ct_45_NLRBC_2nd <- f.roughGate(tp_cp_45_wbc_rough, "UV379-A", thrd = 1.5, tinypeak.removal = 1/25, alpha = 1/40)
  
  
  if (show.plots == T) {
    flowPlot(tp_cp_45_wbc_rough, channels = c("UV379-A", "SSC-A"),
             xlim = c(-1, 4.5), ylim = c(2.5, 5.5), title = identifier(x))
    abline(v = tp_ct_45_NLRBC_2nd)
    
  }
  
  tp_cp_45_NLRBC <- flowDensity(x, c("UV379-A", "SSC-A"), position = c(F, NA), gates = c(tp_ct_45_NLRBC_2nd, NA))
  tp_cp_45_WBC <- flowDensity(x, c("UV379-A", "SSC-A"), position = c(T, NA), gates = c(tp_ct_45_NLRBC_2nd, NA))
  
  if (show.plots == T) {
    flowPlot(x, channels = c("UV379-A", "SSC-A"), xlim = c(-1, 4.5), ylim = c(2.5, 5.5), )
    lines(tp_cp_45_NLRBC@filter, col = "blue", lty = 1, lwd = 2)
    lines(tp_cp_45_WBC@filter, col = "green", lty = 1, lwd = 2)
  }
  
  tp_cp_45_NLRBC@filter[, "SSC-A"] <- 10^tp_cp_45_NLRBC@filter[, "SSC-A"]
  tp_cp_45_WBC@filter[, "SSC-A"] <- 10^tp_cp_45_WBC@filter[, "SSC-A"]
  
  return(list = c(
    NLRBC = polygonGate(filterId = "NLRBC", .gate = tp_cp_45_NLRBC@filter),
    WBC = polygonGate(filterId = "WBC", .gate = tp_cp_45_WBC@filter)
  ))
  
}

















#' @param x a flowrame

#' x = cytoframe_to_flowFrame(fs[[1]]) ; identifier (x); show.plots = T
f.Expanded_Lymph_gating2 <- function(x, # obj,
                                     Tissue, SSC_rot_preSet =1.9,
                                     show.plots = T) {
  # get tissue infor based on the flowframe discription data
  if (missing(Tissue)) {
    Tissue <- c("BM", "PB", "FL", "LN", "OT")

    Tissue <- Tissue[unlist(lapply(Tissue, function(t) {
      grepl(t, identifier(x))
    }))]
  }
  # par(mfrow = c(1, 1))
  # x <-  fs.clean[[1]]
  # cat("performing expanded lymph gate on", identifier(x), "\n")
  # 
  # # 1. FIND BEST NLRBC GATING ON CD45
  # ## 1A. log transform the SSC-A channel for better population identification
   exprs(x)[, "SSC-A"] <- log10(exprs(x)[, "SSC-A"]) # log ssc-A
  # # flowPlot(x, channels = c("UV379-A", "SSC-A"))
  # 
  # # 2. find a cutoff to gate wbc
  # # tried a few ways, the all cut seems works best
  # ct_CD45_UV379_NLRBC_preSet <- 1.5
  # pks <- getPeaks(x, "UV379-A", tinypeak.removal = 1 / 25)
  # 
  # if (sum(pks$Peaks < ct_CD45_UV379_NLRBC_preSet) == 0) {
  #   tp_ct_45_NLRBC <- deGate(x, "UV379-A", use.upper = T, upper = F, tinypeak.removal = 1 / 25, alpha = 1 / 100)
  # } else {
  #   tp_ct_45_NLRBC <- deGate(x, "UV379-A", all.cuts = T, tinypeak.removal = 1 / 25)
  # }
  # 
  # tp_ct_45_NLRBC <- tp_ct_45_NLRBC[which.min(abs(tp_ct_45_NLRBC - ct_CD45_UV379_NLRBC_preSet))]


  if (show.plots == T & nrow(x) >10) {
    par(mfrow = c(1, 4))
    flowPlot(x, channels = c("UV379-A", "SSC-A"), xlim = c(-1, 4.5), ylim = c(2.5, 5.5), title = identifier(x))
    #abline(v = tp_ct_45_NLRBC)
  }


  # Perform gating
  # cat("Set WBC gate at ", tp_ct_45_NLRBC , "\n")
  # tp_cp_45_wbc <- flowDensity(x, c("UV379-A", "SSC-A"), position = c(T, NA), gates = c(tp_ct_45_NLRBC, NA))


  # rotate gate to gate lymph and grans
  rotate.theta <- 0.68
  cat("rotating the data set by", rotate.theta, "\n")

  tp_fs.rot_45SSCA <- rotate.data(x, c("UV379-A", "SSC-A"), theta = rotate.theta)$data
  # flowPlot(tp_fs.rot_45SSCA, c("UV379-A", "SSC-A"), density.overlay = c(T,T)  )

  # Expaned Lymph Gate

  ct_SSC_rot_preSet <- SSC_rot_preSet
  pks <- getPeaks(tp_fs.rot_45SSCA, "SSC-A", tinypeak.removal = 1 / 10)


  tp_ct_SSC_rot <- f.roughGate(tp_fs.rot_45SSCA, "SSC-A", thrd = ct_SSC_rot_preSet)

  #
  # if(sum(pks$Peaks > ct_SSC_rot_preSet) == 0) { # situation that with no Grans
  #   tp_ct_SSC_rot <- deGate(tp_fs.rot_45SSCA, "SSC-A", use.upper = T, upper = T, tinypeak.removal = 1/10)
  # } else {
  #   tp_ct_SSC_rot <- deGate(tp_fs.rot_45SSCA, "SSC-A", all.cuts = T)
  # }
  #
  # tp_ct_SSC_rot <- tp_ct_SSC_rot[which.min(abs(tp_ct_SSC_rot - ct_SSC_rot_preSet))]
  #
  if (show.plots == T & nrow(tp_fs.rot_45SSCA) >10) {
    flowPlot(tp_fs.rot_45SSCA, channels = c("UV379-A", "SSC-A"), title = "first ssc cut") # ,xlim = c(3, 5), ylim = c(2.5, 5.5))
    abline(h = tp_ct_SSC_rot)
  }
  cat("rough gate on rotated WBC at", tp_ct_SSC_rot, "\n")
  # this gate often excluded mono from lymph, so additional gate need to perform on the grans/mono pop
  tp_cp_SSC_grans.rot <- flowDensity(tp_fs.rot_45SSCA, c("UV379-A", "SSC-A"), position = c(NA, T), gates = c(NA, tp_ct_SSC_rot))


  # monos
  if (Tissue %in% c("BM", "PB")) { 
    p <- 1 / 2
    twin.factor = 0.9
  } else {
    p <- 1 / 25
    twin.factor = 0.98
  }
  
  pks <- getPeaks(tp_cp_SSC_grans.rot, "SSC-A", tinypeak.removal = p, twin.factor = twin.factor) # Setting twin factor low to exclude grans in some BM samples

  if (length(pks$Peaks) == 1) { # situation that with no Grans
    tp_ct_SSC_rot <- deGate(tp_cp_SSC_grans.rot, "SSC-A", use.upper = T, upper = F,, twin.factor = 0.6, tinypeak.removal = p) # try to include a bit grans to make sure have all the monos
    tp_ct_SSC_rot <- c(tp_ct_SSC_rot, deGate(tp_cp_SSC_grans.rot, "SSC-A", use.upper = T, upper = T,, twin.factor = 0.6, tinypeak.removal = p) ) # conditions that twin factor merged all peaks
     } else {
    tp_ct_SSC_rot <- deGate(tp_cp_SSC_grans.rot, "SSC-A", all.cuts = T, tinypeak.removal = p)
    # tp_ct_SSC_rot <- c(tp_ct_SSC_rot, deGate(tp_cp_SSC_grans.rot, "SSC-A", use.upper = T, upper = F, tinypeak.removal = 1/25)) #
  }

  tp_ct_SSC_rot <- tp_ct_SSC_rot[which.min(abs(tp_ct_SSC_rot - ct_SSC_rot_preSet))]


  if (show.plots == T & tp_cp_SSC_grans.rot@cell.count >10) {
    flowPlot(tp_cp_SSC_grans.rot, channels = c("UV379-A", "SSC-A")) # ,xlim = c(3, 5), ylim = c(2.5, 5.5))
    abline(h = tp_ct_SSC_rot, col = "red")
  }




 # exp lymph population 
  tp_cp_SSC_expLymph.rot <- flowDensity(tp_fs.rot_45SSCA, c("UV379-A", "SSC-A"), position = c(NA, F), gates = c(NA, tp_ct_SSC_rot))
  # plotDens(tp_cp_SSC_expLymph.rot, c("UV379-A", "SSC-A"), density.overlay = c(T,T)  )

  tp_cp_SSC_expLymph <- rotate.fd(tp_cp_SSC_expLymph.rot, angle = rotate.theta)
  # plotDens(tp_cp_SSC_expLymph, c("UV379-A", "SSC-A"), density.overlay = c(T,T)  )
  tp_cp_SSC_expLymph <- f.remove_pop_outliers(tp_cp_SSC_expLymph, c("UV379-A", "SSC-A"), axis.x = T, axis.y = T)


  # Grans Gate
  tp_cp_SSC_grans.rot <- flowDensity(tp_fs.rot_45SSCA, c("UV379-A", "SSC-A"), position = c(NA, T), gates = c(NA, tp_ct_SSC_rot))
  # plotDens(tp_cp_SSC_expLymph.rot, c("UV379-A", "SSC-A"), density.overlay = c(T,T)  )

  tp_cp_SSC_grans <- rotate.fd(tp_cp_SSC_grans.rot, angle = rotate.theta)

  tp_ct_SSC_grans <- deGate(tp_cp_SSC_grans, "SSC-A", use.upper = T, upper = F)
  tp_cp_SSC_grans <- flowDensity(tp_cp_SSC_grans, c("UV379-A", "SSC-A"), position = c(NA, T), gates = c(NA, tp_ct_SSC_grans))

  

  tp_cp_SSC_expLymph@filter <- f.GatePrettify(tp_cp_SSC_expLymph@filter, ur =T, lr =T)
  tp_cp_SSC_grans@filter <- f.GatePrettify(tp_cp_SSC_grans@filter, ul = T, ur = T )

  if (show.plots == T & nrow(x) > 10) {
    flowPlot(x, channels = c("UV379-A", "SSC-A"), xlim = c(-1, 4.5), ylim = c(3.5, 5.5), )
    lines(tp_cp_SSC_expLymph@filter, col = "blue", lty = 1, lwd = 2)
    lines(tp_cp_SSC_grans@filter, col = "green", lty = 1, lwd = 2)
  }

  exprs(x)[, "SSC-A"] <- 10^(exprs(x)[, "SSC-A"]) # log ssc-A
  tp_cp_SSC_expLymph@filter[, "SSC-A"] <- 10^tp_cp_SSC_expLymph@filter[, "SSC-A"]
  tp_cp_SSC_grans@filter[, "SSC-A"] <- 10^tp_cp_SSC_grans@filter[, "SSC-A"]
  # Translate the filter to gate



  return(list = c(
    Grans = polygonGate(filterId = "Grans", .gate = tp_cp_SSC_grans@filter),
    Exp_lymph = polygonGate(filterId = "Exp_Lymph", .gate = tp_cp_SSC_expLymph@filter)
  ))
}



# f.roughGate <- function(x, channel, thrd.x, thrd.y,  tinypeak.removal = 1/50, afterpeak  = T, alpha = 1/10) {
#
#   pks <- getPeaks(x, channel, tinypeak.removal = tinypeak.removal)
#
#
#   if(missing(thrd.x)) { peak.condition = (length(pks$Peaks) ==1) } else {
#     peak.condition  <- (sum(pks$Peak > thrd.x) ==0 )
#   }
#
#   if(peak.condition) {
#     cut <- deGate(x, channel, use.upper = T, upper = afterpeak, tinypeak.removal = tinypeak.removal, alpha = alpha)
#   } else {
#     cut <- deGate(x, channel, tinypeak.removal = tinypeak.removal)
#   }
#   return(cut)
# }











#' pre PDC_baso
#' x =  cytoframe_to_flowFrame(fs[[1]]); identifier(x); show.plots = T; x <- cytoframe_to_flowFrame(x)
f.prePDC_Baso <- function(x,
                       show.plots = T) {
  
  
  sample_name <- identifier(x)
  print(sample_name)
  
  
  cat("STEP 1: rough gating on CD19(B750) to remove potential CD123 B cells(hairy cells)\n")
  tp_ct_B750 <- f.roughGate(x, "B750-A", tinypeak.removal = 1 / 100, thrd = 2.5, alpha = 1 / 2)
  if (show.plots == T & nrow(x) >10) {
    par(mfrow = c(1, 4))
    flowPlot(x, channels = c("B750-A", "YG730-A"), title = identifier(x), xlim = c(0, 4.5), ylim = c(0, 4.5))
    abline( v = tp_ct_B750)
  }
  

  tp_cp_YG730_T <- flowDensity(x, c("B750-A", "YG730-A"), position = c(F, NA), gates = c(tp_ct_B750, NA))

  
  
  
  l.CD19n <- f.generate_gate(list(`CD19nWBC` = tp_cp_YG730_T), gate.type = "r")

  
  l.CD19n_cood <- f.extractGate_boundary(l.CD19n)
  
  
  
  if (show.plots == T & tp_cp_YG730_T@cell.count >10) {
    flowPlot(x, channels = c("B750-A", "YG730-A"), title = "refine 660 cut", xlim = c(0, 4.5), ylim = c(0, 4.5))
    lapply(l.CD19n_cood, function(x) {
      lines(x, col = "red", lty = 1, lwd = 2)
    })
  }
  
  
  return(l.CD19n)

}









#' PDC Baso
#' x =  cytoframe_to_flowFrame(fs[[1]]); identifier(x); show.plots = T; x <- cytoframe_to_flowFrame(x)
f.PDC_Baso <- function(x,
                       show.plots = T) {
  
  
  sample_name <- identifier(x)
  print(sample_name)
  
  # if(grepl("QC", sample_name) & grepl("CHEX", sample_name) & grepl("103", sample_name)) {
  #   warning("No PDC/Baso gate for CD-Chex CD103 Plus QC Sample")
  #   
  #   return(NULL)
  # } else {
    cat("STEP 1: rough gating on CD123(YG730) to remove main peak\n")
    tp_ct_YG730 <- f.roughGate(x, "YG730-A", tinypeak.removal = 1 / 10, thrd = 2.5, alpha = 1 / 50)
    #tp_ct_SSC <- f.roughGate(x, "V510-A", tinypeak.removal = 1 / 100, thrd = 4.5)
    
    # remove main peaks
    tp_cp_YG730_T <- flowDensity(x, c("YG730-A", "V510-A"), position = c(T, NA), gates = c(tp_ct_YG730, NA))
    
    cat("STEP 2: finding gate for PDC and BASO\n")
    
    tp_ct_YG730_2 <- f.roughGate(tp_cp_YG730_T, "YG730-A", tinypeak.removal = 1 / 100, thrd = 3, alpha = 1 / 50)
    
    tp_cp_YG730 <- flowDensity(x, c("YG730-A", "V510-A"), position = c(T, NA), gates = c(tp_ct_YG730_2, NA))
    
    
    if (show.plots == T & tp_cp_YG730_T@cell.count >10) {
      par(mfrow = c(1, 4))
      flowPlot(tp_cp_YG730_T, channels = c("YG730-A", "V510-A"), title = identifier(x), xlim = c(0, 4.5), ylim = c(0, 4.5))
      abline( v = tp_ct_YG730)
      abline( v = tp_ct_YG730_2)
      # abline(v = tp_ct_YG660_lymph, col = "red")
    }
    # flowPlot(tp_cp_14SSCA_lymph_pre, channels = c("YG660-A", "SSC-A"))
    # flowPlot(tp_cp_14SSCA_mono_pre, channels = c("YG660-A", "SSC-A"))
    
    cat("STEP 3: rough V510 gate\n")
    # second round of gate
    #tp_ct_YG660_lymph <- f.roughGate(tp_cp_14SSCA_lymph_pre, "YG660-A", tinypeak.removal = 1 / 50, thrd = 2, alpha = 1 / 20) # incread alpha to make cut tighter for the second round
    tp_ct_DR <- f.roughGate(tp_cp_YG730, "V510-A", tinypeak.removal = 1 / 20, thrd = 1.5, alpha = 1/50) # push a bit
    
    
    if (show.plots == T & tp_cp_YG730@cell.count >10) {
      flowPlot(tp_cp_YG730, channels = c("YG730-A", "V510-A"), title = "refine 660 cut", xlim = c(0, 4.5), ylim = c(0, 4.5))
      abline(h = tp_ct_DR, v = tp_ct_YG730_2, col = "red")
      # abline(v = tp_ct_YG660_lymph, col = "red")
    }
    
    
    
    tp_cp_PDC <- flowDensity(tp_cp_YG730, c("YG730-A", "V510-A"), position = c(NA, T), gates = c(NA, tp_ct_DR)) 
    
    
    cat("STEP 4: Fine tuning V510 gate\n")
    tp_ct_PDC <- max(tp_ct_DR, f.roughGate(tp_cp_PDC, "V510-A", tinypeak.removal = 1 / 100, thrd = 1.5, alpha = 1/10)) # push a bit
    
    # get baso 
    tp_cp_Baso <- flowDensity(tp_cp_YG730_T, c("YG730-A", "V510-A"), position = c(NA, F), gates = c(NA, tp_ct_PDC)) 
    cat("Basp\n")
    tp_ct_baso <- min(tp_ct_DR, f.roughGate(tp_cp_Baso, "V510-A", tinypeak.removal = 1 / 100, thrd = 1.5, alpha = 1/10)) # push a bit
    cat(tp_ct_baso)
    tp_ct_123 <- f.roughGate(tp_cp_Baso, "YG730-A", tinypeak.removal = 1 / 100, thrd = 3, alpha = 1/10)
    
    PDC <- flowDensity(tp_cp_YG730, c("YG730-A", "V510-A"), position = c(NA, T), gates = c(NA, tp_ct_PDC)) 
    Baso <- flowDensity(tp_cp_Baso, c("YG730-A", "V510-A"), position = c(T, F), gates = c(tp_ct_123, tp_ct_baso)) 
    
    
    
    
    l.pdcbaso <- f.generate_gate(list(`PDC` = PDC,
                                      `Baso` = Baso
    ), gate.type = "r")
    # l.56and3_gate$NK@max[2] <- 4.3
    l.pdcbaso_gate_cood <- f.extractGate_boundary(l.pdcbaso)
    
    
    
    if (show.plots == T & tp_cp_Baso@cell.count >10) {
      flowPlot(tp_cp_Baso, channels = c("YG730-A", "V510-A"), title = "refine 660 cut", xlim = c(0, 4.5), ylim = c(0, 4.5))
      abline(h = tp_ct_PDC, v = tp_ct_YG730, col = "green")
      abline(h = tp_ct_baso, v = tp_ct_YG730, col = "blue")
      # abline(v = tp_ct_YG660_lymph, col = "red")
    }
    
    if (show.plots == T & nrow(x) >10) {
      flowPlot(x, channels = c("YG730-A", "V510-A"), title = "refine 660 cut", xlim = c(0, 4.5), ylim = c(0, 4.5))
      lapply(l.pdcbaso_gate_cood, function(x) {
        lines(x, col = "red", lty = 1, lwd = 2)
      })
      
      
    }
    
    
    return(l.pdcbaso)
  # }
  

  
    # flowPlot(x, channels = c("YG660-A", "SSC-A"))
  

  # return(list = c(
  #   Lymph = polygonGate(filterId = "Lymph", .gate = tp_cp_14SSCA_lymph@filter),
  #   Mono = polygonGate(filterId = "Mono", .gate = tp_cp_14SSCA_mono@filter)
  # ))
}




















#' Mono gating
#' x = fs[[1]]; identifier(x); show.plots = T; x <- cytoframe_to_flowFrame(x)
f.lymph_mono_gating_B <- function(x,
                                  show.plots = T) {
  print(identifier(x))
  exprs(x)[, "SSC-A"] <- log10(exprs(x)[, "SSC-A"]) # log ssc-A
  # flowPlot(x, channels = c("YG660-A", "SSC-A"))
  cat("STEP 1: rough gating\n")
  tp_ct_YG660 <- f.roughGate(x, "YG660-A", tinypeak.removal = 1 / 50, thrd = 2, alpha = 1 / 50)
  tp_ct_SSC <- f.roughGate(x, "SSC-A", tinypeak.removal = 1 / 50, thrd = 4.5)

  if (show.plots == T & nrow(x) >5) {
    par(mfrow = c(1, 4))
    flowPlot(x, channels = c("YG660-A", "SSC-A"), title = identifier(x), xlim = c(0, 4.5), ylim = c(3.5, 5.5))
    abline(h = tp_ct_SSC, v = tp_ct_YG660)
    # abline(v = tp_ct_YG660_lymph, col = "red")
  }

  tp_cp_14SSCA_lymph_pre <- flowDensity(x, c("YG660-A", "SSC-A"), position = c(F, NA), gates = c(tp_ct_YG660, NA))
  tp_cp_14SSCA_mono_pre <- flowDensity(x, c("YG660-A", "SSC-A"), position = c(T, NA), gates = c(tp_ct_YG660, NA))
  # flowPlot(tp_cp_14SSCA_lymph_pre, channels = c("YG660-A", "SSC-A"))
  # flowPlot(tp_cp_14SSCA_mono_pre, channels = c("YG660-A", "SSC-A"))





 cat("STEP 2: refine SSC gate\n")
  # second round of gate
  #tp_ct_YG660_lymph <- f.roughGate(tp_cp_14SSCA_lymph_pre, "YG660-A", tinypeak.removal = 1 / 50, thrd = 2, alpha = 1 / 20) # incread alpha to make cut tighter for the second round
  tp_ct_SSC_lymph <- f.roughGate(tp_cp_14SSCA_lymph_pre, "SSC-A", tinypeak.removal = 1 / 10, thrd = 5, alpha = 1/50) # push a bit
  tp_ct_SSC_mono <- f.roughGate(tp_cp_14SSCA_mono_pre, "SSC-A", tinypeak.removal = 1 / 50, thrd = 4.3)
  tp_ct_YG660_mono <- max(f.roughGate(tp_cp_14SSCA_mono_pre, "YG660-A", tinypeak.removal = 1 / 50, alpha = 1 / 4), tp_ct_YG660)

  if (show.plots == T & tp_cp_14SSCA_lymph_pre@cell.count >5) {
    flowPlot(tp_cp_14SSCA_lymph_pre, channels = c("YG660-A", "SSC-A"), title = "refine 660 cut", xlim = c(0, 4.5), ylim = c(3.5, 5.5))
    abline(h = tp_ct_SSC_lymph, v = tp_ct_YG660, col = "red")
    # abline(v = tp_ct_YG660_lymph, col = "red")
  }


  
  
  tp_cp_14SSCA_lymph <- flowDensity(x, c("YG660-A", "SSC-A"), position = c(F, F), gates = c(tp_ct_YG660, tp_ct_SSC_lymph)) # 250210: tp_ct_YG660_lymph changed to tp_ct_YG660 to accomondate potential undercomp issue
  
  if(is.infinite(tp_ct_SSC_mono) | is.infinite(tp_ct_YG660_mono)) {
    tp_cp_14SSCA_mono <- tp_cp_14SSCA_mono_pre
  } else {
    tp_cp_14SSCA_mono <- flowDensity(x, c("YG660-A", "SSC-A"), position = c(T, T), gates = c(tp_ct_YG660_mono, tp_ct_SSC_mono))
    
  }
  
  
  if (show.plots == T & tp_cp_14SSCA_mono_pre@cell.count >5) {
    flowPlot(tp_cp_14SSCA_mono_pre, channels = c("YG660-A", "SSC-A"), title = identifier(x), xlim = c(0, 4.5), ylim = c(3.5, 5.5))
    abline(h = tp_ct_SSC_mono, v = tp_ct_YG660_mono, col = "green")
    # abline(v = tp_ct_YG660_lymph, col = "red")
  }

  # outlier_cut_x<-deGate(tp_cp_14SSCA_lymph,"YG660-A", use.percentile = T, percentile = 0.0001)
  # outlier_cut_y<-deGate(tp_cp_14SSCA_lymph,"SSC-A", use.percentile = T, percentile = 0.0001)


  tp_cp_14SSCA_lymph <- f.remove_pop_outliers(tp_cp_14SSCA_lymph, c("YG660-A", "SSC-A"), axis.x = T, axis.y = T)
  tp_cp_14SSCA_mono <- f.remove_pop_outliers(tp_cp_14SSCA_mono, c("YG660-A", "SSC-A"), axis.x = T, axis.y = T)

  tp_cp_14SSCA_lymph@filter <- f.GatePrettify(tp_cp_14SSCA_lymph@filter, ll = T, ul = T)
  if(tp_cp_14SSCA_mono@cell.count >20) {
  tp_cp_14SSCA_mono@filter <- f.GatePrettify(tp_cp_14SSCA_mono@filter, ur = T, lr = T, ul =T)
  }

  # tp_cp_14SSCA_lymph <- flowDensity(tp_cp_14SSCA_lymph, c("YG660-A", "SSC-A"), position = c(T, T), gates = c(outlier_cut_x, outlier_cut_y))


  if (show.plots == T & nrow(x) >5) {
    flowPlot(x, channels = c("YG660-A", "SSC-A"), xlim = c(0, 4.5), ylim = c(3.5, 5.5))
    lines(tp_cp_14SSCA_lymph@filter, col = "blue", lty = 1, lwd = 2)
    lines(tp_cp_14SSCA_mono@filter, col = "green", lty = 1, lwd = 2)
  }

  exprs(x)[, "SSC-A"] <- 10^(exprs(x)[, "SSC-A"]) # log ssc-A
  tp_cp_14SSCA_lymph@filter[, "SSC-A"] <- 10^tp_cp_14SSCA_lymph@filter[, "SSC-A"]
  tp_cp_14SSCA_mono@filter[, "SSC-A"] <- 10^tp_cp_14SSCA_mono@filter[, "SSC-A"]


  return(list = c(
    Mono = polygonGate(filterId = "Mono", .gate = tp_cp_14SSCA_mono@filter),
    Lymph = polygonGate(filterId = "Lymph", .gate = tp_cp_14SSCA_lymph@filter)

  ))
}





#' Mono gating
#' x = fs[[1]]; identifier(x); show.plots = T
f.lymph_mono_gating_T <- function(x,ssc_preset = 4.5,
                                  show.plots = T) {
  print(identifier(x))
  
  if(nrow(x) <5) {
    
    dummy_gate <- matrix(rep(0,6), nrow = 3)
    colnames(dummy_gate) <- c("B602-A", "SSC-A")
    
    return(list = c(
      Lymph = polygonGate(filterId = "Lymph", .gate = dummy_gate),
      Mono = polygonGate(filterId = "Mono", .gate = dummy_gate)
    ))
  }else {
    
    exprs(x)[, "SSC-A"] <- log10(exprs(x)[, "SSC-A"]) # log ssc-A
    # flowPlot(x, channels = c("B602-A", "SSC-A"))
    tp_ct_B602 <- f.roughGate(x, "B602-A", tinypeak.removal = 1 / 50, thrd = 2, alpha = 1 / 50)
    tp_ct_SSC <- f.roughGate(x, "SSC-A", tinypeak.removal = 1 / 50, thrd = ssc_preset)
    
    if (show.plots == T) {
      par(mfrow = c(1, 4))
      flowPlot(x, channels = c("B602-A", "SSC-A"), title = identifier(x), xlim = c(0, 4.5), ylim = c(3.5, 5.5))
      abline(h = tp_ct_SSC, v = tp_ct_B602)
      # abline(v = tp_ct_B602_lymph, col = "red")
    }
    
    tp_cp_14SSCA_lymph_pre <- flowDensity(x, c("B602-A", "SSC-A"), position = c(F, NA), gates = c(tp_ct_B602, NA))
    tp_cp_14SSCA_mono_pre <- flowDensity(x, c("B602-A", "SSC-A"), position = c(T, NA), gates = c(tp_ct_B602, NA))
    # flowPlot(tp_cp_14SSCA_lymph_pre, channels = c("B602-A", "SSC-A"))
    # flowPlot(tp_cp_14SSCA_mono_pre, channels = c("B602-A", "SSC-A"))
  
    
    # second round of gate
    tp_ct_B602_lymph <- f.roughGate(tp_cp_14SSCA_lymph_pre, "B602-A", tinypeak.removal = 1 / 50, thrd = 2, alpha = 1 / 20) # incread alpha to make cut tighter for the second round
    tp_ct_SSC_lymph <- f.roughGate(tp_cp_14SSCA_lymph_pre, "SSC-A", tinypeak.removal = 1 / 10, thrd = ssc_preset)
    
    
    
    
    tp_ct_SSC_mono <- f.roughGate(tp_cp_14SSCA_mono_pre, "SSC-A", tinypeak.removal = 1 / 50, thrd = ssc_preset-0.5)
  
    tp_ct_B602_mono <- f.roughGate(tp_cp_14SSCA_mono_pre, "B602-A", tinypeak.removal = 1 / 50,  thrd = 2, alpha = 1 / 4)
    
    
    if (show.plots == T & tp_cp_14SSCA_lymph_pre@cell.count >10) {
      flowPlot(tp_cp_14SSCA_lymph_pre, channels = c("B602-A", "SSC-A"), title = identifier(x), xlim = c(0, 4.5), ylim = c(3.5, 5.5))
      abline(h = tp_ct_SSC_lymph, v = tp_ct_B602_lymph, col = "red")
      # abline(v = tp_ct_B602_lymph, col = "red")
    }
    
    tp_ct_SSC_mono <- min(tp_ct_SSC_mono, tp_ct_SSC)
    tp_ct_B602_mono <- max(tp_ct_B602_mono, tp_ct_B602)
    
    tp_ct_B602_lymph <- min(tp_ct_B602_lymph, tp_ct_B602)
    tp_ct_SSC_lymph <- min(tp_ct_SSC_lymph, tp_ct_SSC)
    
    tp_cp_14SSCA_lymph <- flowDensity(x, c("B602-A", "SSC-A"), position = c(F, F), gates = c(tp_ct_B602_lymph, tp_ct_SSC_lymph))
    tp_cp_14SSCA_mono <- flowDensity(x, c("B602-A", "SSC-A"), position = c(T, T), gates = c(tp_ct_B602_mono, tp_ct_SSC_mono))
    
    
    

    
    if (show.plots == T & tp_cp_14SSCA_mono_pre@cell.count >10) { 
      flowPlot(tp_cp_14SSCA_mono_pre, channels = c("B602-A", "SSC-A"), title = identifier(x), xlim = c(0, 4.5), ylim = c(3.5, 5.5))
      abline(h = tp_ct_SSC_mono, v = tp_ct_B602_mono, col = "green")
      # abline(v = tp_ct_B602_lymph, col = "red")
    }
    
    
    # outlier_cut_x<-deGate(tp_cp_14SSCA_lymph,"B602-A", use.percentile = T, percentile = 0.0001)
    # outlier_cut_y<-deGate(tp_cp_14SSCA_lymph,"SSC-A", use.percentile = T, percentile = 0.0001)
    
    
    tp_cp_14SSCA_lymph <- f.remove_pop_outliers(tp_cp_14SSCA_lymph, c("B602-A", "SSC-A"), axis.x = T, axis.y = T)
    
    tp_cp_14SSCA_mono <- f.remove_pop_outliers(tp_cp_14SSCA_mono, c("B602-A", "SSC-A"), axis.x = T, axis.y = T)
    
  
    
    tp_cp_14SSCA_lymph@filter <- f.GatePrettify(tp_cp_14SSCA_lymph@filter, ll = T, ul = T)
    
    if(tp_cp_14SSCA_mono@cell.count >20) {
      tp_cp_14SSCA_mono@filter <- f.GatePrettify(tp_cp_14SSCA_mono@filter, ur = T, lr = T, ul =T)
    }
    
    # tp_cp_14SSCA_lymph <- flowDensity(tp_cp_14SSCA_lymph, c("B602-A", "SSC-A"), position = c(T, T), gates = c(outlier_cut_x, outlier_cut_y))

    if (show.plots == T & nrow(x) >5) {
      flowPlot(x, channels = c("B602-A", "SSC-A"), xlim = c(0, 4.5), ylim = c(3.5, 5.5))
      lines(tp_cp_14SSCA_lymph@filter, col = "blue", lty = 1, lwd = 2)
      lines(tp_cp_14SSCA_mono@filter, col = "green", lty = 1, lwd = 2)
    }
    
    exprs(x)[, "SSC-A"] <- 10^(exprs(x)[, "SSC-A"]) # log ssc-A
    tp_cp_14SSCA_lymph@filter[, "SSC-A"] <- 10^tp_cp_14SSCA_lymph@filter[, "SSC-A"]
    tp_cp_14SSCA_mono@filter[, "SSC-A"] <- 10^tp_cp_14SSCA_mono@filter[, "SSC-A"]
    
    
    return(list = c(
      Mono = polygonGate(filterId = "Mono", .gate = tp_cp_14SSCA_mono@filter),
      Lymph = polygonGate(filterId = "Lymph", .gate = tp_cp_14SSCA_lymph@filter)

    ))
    
  }
  
  
  
  
  
}




#' Mono gating
#' x = gs_pop_get_data(gs, "Exp_Lymph")[[1]]
f.lymph_mono_gating <- function(x,
                                show.plots = T) {
  print(identifier(x))
  exprs(x)[, "SSC-A"] <- log10(exprs(x)[, "SSC-A"]) # log ssc-A
  # flowPlot(x, channels = c("B602-A", "SSC-A"))


  # lots of negative events, so perform a clean up gate
  tp_ct_B602 <- deGate(x, "B602-A", use.upper = T, upper = F)
  tp_cp_B602_clean <- flowDensity(x, c("B602-A", "SSC-A"), position = c(T, NA), gates = c(tp_ct_B602, NA))



  rotate.theta <- -0.5
  tp_fs.rot_14SSCA <- rotate.data(tp_cp_B602_clean@flow.frame, c("B602-A", "SSC-A"), theta = rotate.theta)$data
  # plotDens(tp_fs.rot_14SSCA, c("B602-A", "SSC-A"), density.overlay = c(T,T)  )

  tp_pk_SSC.rot <- getPeaks(tp_fs.rot_14SSCA, "SSC-A")

  if (length(tp_pk_SSC.rot$Peaks) == 1) {
    tp_ct_SSC.rot <- deGate(tp_fs.rot_14SSCA, "SSC-A", upper = T)
  } else {
    tp_ct_SSC.rot <- deGate(tp_fs.rot_14SSCA, "SSC-A")
  }

  tp_cp_14SSCA_lymph.rot <- flowDensity(tp_fs.rot_14SSCA, c("B602-A", "SSC-A"), position = c(NA, F), gates = c(NA, tp_ct_SSC.rot))
  tp_cp_14SSCA_mono.rot <- flowDensity(tp_fs.rot_14SSCA, c("B602-A", "SSC-A"), position = c(NA, T), gates = c(NA, tp_ct_SSC.rot))

  tp_cp_14SSCA_lymph <- rotate.fd(tp_cp_14SSCA_lymph.rot, angle = rotate.theta)
  tp_cp_14SSCA_mono <- rotate.fd(tp_cp_14SSCA_mono.rot, angle = rotate.theta)

  if (show.plots == T) {
    par(mfrow = c(1, 4))
    flowPlot(x, channels = c("B602-A", "SSC-A"), title = identifier(x), xlim = c(0, 4.5))
    plotDens(tp_fs.rot_14SSCA, channels = c("B602-A", "SSC-A"))
    abline(h = tp_ct_SSC.rot)
    flowPlot(x, channels = c("B602-A", "SSC-A"), xlim = c(0, 4.5))
    lines(tp_cp_14SSCA_lymph@filter, col = "blue", lty = 1, lwd = 2)
    lines(tp_cp_14SSCA_mono@filter, col = "green", lty = 1, lwd = 2)
  }

  exprs(x)[, "SSC-A"] <- 10^(exprs(x)[, "SSC-A"]) # log ssc-A
  tp_cp_14SSCA_lymph@filter[, "SSC-A"] <- 10^tp_cp_14SSCA_lymph@filter[, "SSC-A"]
  tp_cp_14SSCA_mono@filter[, "SSC-A"] <- 10^tp_cp_14SSCA_mono@filter[, "SSC-A"]


  return(list = c(
    Lymph = polygonGate(filterId = "Lymph", .gate = tp_cp_14SSCA_lymph@filter),
    Mono = polygonGate(filterId = "Mono", .gate = tp_cp_14SSCA_mono@filter)
  ))
}



#' TB gating
#' x <- gs_pop_get_data(gs, "Lymph")[[3]]


f.T_B_gating <- function(x,
                         show.plots = T) {
  rotate.theta <- -0.9
  tp_fs.rot_CD3CD19 <- rotate.data(x, c("B810-A", "B750-A"), theta = rotate.theta)$data
  # plotDens(tp_fs.rot_CD3CD19, c("B810-A", "B750-A"), density.overlay = c(T,T)  )

  tp_ct_19.rot <- deGate(tp_fs.rot_CD3CD19, "B750-A")


  tp_cp_CD3CD19_lymph.rot <- flowDensity(tp_fs.rot_CD3CD19, c("B810-A", "B750-A"), position = c(NA, F), gates = c(NA, tp_ct_19.rot))
  tp_cp_CD3CD19_TandB.rot <- flowDensity(tp_fs.rot_CD3CD19, c("B810-A", "B750-A"), position = c(NA, T), gates = c(NA, tp_ct_19.rot))

  # plotDens(tp_cp_CD3CD19_TandB.rot, c("B810-A", "B750-A"), density.overlay = c(T,T)  )

  tp_ct_3.rot <- deGate(tp_cp_CD3CD19_TandB.rot, "B810-A")


  tp_cp_CD3CD19_B.rot <- flowDensity(tp_cp_CD3CD19_TandB.rot, c("B810-A", "B750-A"), position = c(F, NA), gates = c(tp_ct_3.rot, NA))
  tp_cp_CD3CD19_T.rot <- flowDensity(tp_cp_CD3CD19_TandB.rot, c("B810-A", "B750-A"), position = c(T, NA), gates = c(tp_ct_3.rot, NA))

  tp_cp_CD3CD19_lymph <- rotate.fd(tp_cp_CD3CD19_lymph.rot, angle = rotate.theta)
  tp_cp_CD3CD19_B <- rotate.fd(tp_cp_CD3CD19_B.rot, angle = rotate.theta)
  tp_cp_CD3CD19_T <- rotate.fd(tp_cp_CD3CD19_T.rot, angle = rotate.theta)

  # refine lymph gate by cutting the CD19 "pos"
  tp_ct_19 <- deGate(x, "B750-A")
  tp_cp_CD3CD19_lymph <- flowDensity(tp_cp_CD3CD19_lymph, c("B810-A", "B750-A"), position = c(NA, F), gates = c(NA, tp_ct_19))


  if (show.plots == T) {
    par(mfrow = c(1, 4))
    flowPlot(x, channels = c("B810-A", "B750-A"), xlim = c(0, 4.5))
    flowPlot(tp_fs.rot_CD3CD19, channels = c("B810-A", "B750-A"))
    flowPlot(x, channels = c("B810-A", "B750-A"), xlim = c(0, 4.5))
    lines(tp_cp_CD3CD19_lymph@filter, col = "blue", lty = 1, lwd = 2)
    lines(tp_cp_CD3CD19_B@filter, col = "green", lty = 1, lwd = 2)
    lines(tp_cp_CD3CD19_T@filter, col = "cyan", lty = 1, lwd = 2)
  }

  return(list = c(
    B_cell = polygonGate(filterId = "B_cell", .gate = tp_cp_CD3CD19_B@filter),
    T_cell = polygonGate(filterId = "T_cell", .gate = tp_cp_CD3CD19_T@filter),
    NonBT = polygonGate(filterId = "NonBT", .gate = tp_cp_CD3CD19_lymph@filter)
  ))
}


#' TB gating
#' x <- gs_pop_get_data(gs, "CD3pCD2p")[[3]]


f.gdT_gating <- function(x,
                         show.plots = T) {
  # rotate.theta <- -0.6
  # tp_fs.rot_CD3TCRgd <- rotate.data(x, c("B810-A", "UV660-A"), theta = rotate.theta)$data
  # plotDens(tp_fs.rot_CD3TCRgd, c("B810-A", "UV660-A"), density.overlay = c(T,T)  )

  # getPeaks(x, "UV660-A", tinypeak.removal = 1/100)

  tp_ct_UV660 <- deGate(x, "UV660-A", tinypeak.removal = 1 / 100, all.cuts = T)

  dens_UV660 <- f.dens_estimation(x, "UV660-A")

  dens.tp_ct_UV660 <- f.get_cutoff_dens(dens_UV660, tp_ct_UV660)

  tp_ct_UV660 <- tp_ct_UV660[which.min(dens.tp_ct_UV660)]

  tp_cp_CD3Tab.rot <- flowDensity(x, c("B810-A", "UV660-A"), position = c(NA, F), gates = c(NA, tp_ct_UV660.rot))
  tp_cp_CD3Tgd.rot <- flowDensity(x, c("B810-A", "UV660-A"), position = c(NA, T), gates = c(NA, tp_ct_UV660.rot))



  if (show.plots == T) {
    par(mfrow = c(1, 4))
    flowPlot(x, channels = c("B810-A", "UV660-A"))
    abline(h = tp_ct_UV660, col = "blue")
  }

  return(list = c(
    gdT = polygonGate(filterId = "gdT_cell", .gate = tp_cp_CD3Tab@filter),
    abT = polygonGate(filterId = "abT_cell", .gate = tp_cp_CD3Tgd@filter),
  ))
}






# x <- gs_pop_get_data(gs, "NonB")[[1]] %>% cytoframe_to_flowFrame()
f.gdT_gating_rot <- function(x,
                             show.plots = T) {
  
  negative_noise_cut <- min(deGate(x, "UV660-A", use.percentile = T, percentile = 0.001), 0)
  x <- getflowFrame(flowDensity(x, c("B810-A", "UV660-A"), position = c(NA, T), gates = c(NA, negative_noise_cut)))
  
  
  # gate to remove the CD3 neg populaiton as it can affect set the gating on UV660 

  
  

  rotate.theta <- -0.2 #-0.6
  tp_fs.rot_CD3TCRgd.rot <- rotate.data(x, c("B810-A", "UV660-A"), theta = rotate.theta)$data
  # flowPlot(tp_fs.rot_CD3TCRgd.rot , c("B810-A", "UV660-A"), density.overlay = c(T, T))

  tp_ct_B810.rot<- f.roughGate(tp_fs.rot_CD3TCRgd.rot, "B810-A", thrd = 1.5, tinypeak.removal = 1 / 1000)
  cat("Find B810-A gate at ", tp_ct_B810.rot, "\n")
  
  tp_fs.rot_CD3TCRgd.rot.B810p <- flowDensity(tp_fs.rot_CD3TCRgd.rot, c("B810-A", "UV660-A"), position = c(T, NA), gates = c(tp_ct_B810.rot, NA))
  
  
  
  tp_ct_UV660.rot<- f.roughGate(tp_fs.rot_CD3TCRgd.rot.B810p, "UV660-A", thrd = 2.5, tinypeak.removal = 1 / 1000)
  cat("Find UV660-A gate at ", tp_ct_UV660.rot, "\n")
  
  
  
  
  
  tp_ct_UV660.rot <- tp_ct_UV660.rot[which.min(abs(tp_ct_UV660.rot - 2.5))]

  tp_cp_CD3Tab.rot <- flowDensity(tp_fs.rot_CD3TCRgd.rot, c("B810-A", "UV660-A"), position = c(NA, F), gates = c(NA, tp_ct_UV660.rot))
  tp_cp_CD3Tgd.rot <- flowDensity(tp_fs.rot_CD3TCRgd.rot, c("B810-A", "UV660-A"), position = c(T, T), gates = c(tp_ct_B810.rot, tp_ct_UV660.rot))


  tp_cp_CD3Tab <- rotate.fd(tp_cp_CD3Tab.rot, angle = rotate.theta)
  tp_cp_CD3Tgd <- rotate.fd(tp_cp_CD3Tgd.rot, angle = rotate.theta)
  
  tp_cp_CD3Tab@filter <- f.GatePrettify(tp_cp_CD3Tab@filter, ll = T, lr = T, ul =T)
  tp_cp_CD3Tgd@filter <- f.GatePrettify(tp_cp_CD3Tgd@filter, ul = T, ur = T, lr =T)
  

  if (show.plots == T) {
    par(mfrow = c(1, 4))
    flowPlot(x, channels = c("B810-A", "UV660-A"), title = identifier(x), xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5))
    flowPlot(tp_fs.rot_CD3TCRgd.rot, channels = c("B810-A", "UV660-A"), xlim = c(0.5, 4.5), ylim = c(1, 4))
    abline(h = tp_ct_UV660.rot, v =tp_ct_B810.rot)
    flowPlot(x, channels = c("B810-A", "UV660-A"), xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5))
    lines(tp_cp_CD3Tab@filter, col = "blue", lty = 1, lwd = 2)
    lines(tp_cp_CD3Tgd@filter, col = "green", lty = 1, lwd = 2)
  }



  return(list = c(
    gdT = polygonGate(filterId = "gdT", .gate = tp_cp_CD3Tgd@filter),
    abT = polygonGate(filterId = "NonGD", .gate = tp_cp_CD3Tab@filter)
  ))
}


# a <- lapply( seq(-1,0, 0.1), function(theta) {
#
#   rot  <- rotate.data(x, c("B810-A", "UV660-A"), theta = theta)
#   plotDens(rot$data , c("B810-A", "UV660-A"), density.overlay = c(T, T), main = theta)
#
#   cutoffs.rot <- deGate(rot$data, "UV660-A", tinypeak.removal = 1 / 100, all.cuts = T)
#
#   dens.rot <- f.dens_estimation(rot$data, "UV660-A")
#
#   dens.cutoffs <- f.get_cutoff_dens(dens.rot, cutoffs.rot)
#
#   lowest_cutoff <- cutoffs.rot[which.min(dens.cutoffs)]
#
#   return(data.frame(theta = theta, dens = min(dens.cutoffs), cutoff =lowest_cutoff) )
#
# }
# )
#
# do.call(rbind,a)


#' @param obj a flowframe
#' @param channel, a marker


f.dens_estimation <- function(obj, channel, adjust.dens = 1, spar = .4) {
  data <- exprs(obj)[, channel]



  if (length(data) < 2) {
    warning("Less than 2 cells, returning NA.")
    return(NA)
  }
  dens <- density(data[which(!is.na(data))], adjust = adjust.dens)
  dens <- smooth.spline(x = dens$x, y = dens$y, spar = spar)
  dens$y[which(dens$y < 0)] <- 0

  return(dens)
}

f.get_cutoff_dens <- function(dens, cutoffs) {
  unlist(lapply(cutoffs, function(x) {
    dens$y[which.min(abs(dens$x - x))] * 100
  }))
}



#' select B cells and NonB cell
#' NonB cells then can be subject to gdT cells then  CD2 and CD3 gating
#'  x<- gs_pop_get_data(gs, "Lymph")[[3]]
#'
f.NonB_gating <- function(x, gate.type = "r",
                          show.plots = T) {
  # plotDens(x, c("B750-A", "SSC-A"), density.overlay = c(T,T)  )

  pks <- getPeaks(x, "B750-A")

  if (length(pks$Peaks) == 1) {
    tp_ct_19 <- deGate(x, "B750-A", use.upper = T, upper = T)
  } else {
    tp_ct_19 <- deGate(x, "B750-A", all.cuts = T)
  }


  f.nearest_to_reference <- function(cuts, reference) {
    cuts[which.min(abs(cuts - reference))]
  }
  
  tp_ct_19 <- f.nearest_to_reference(tp_ct_19, 1.8)
  
  if (tp_ct_19 > 3 |tp_ct_19 < 1 ) { # when there are a lot of B cells
    tp_ct_19 <- 1.8
  }


  tp_cp_19_B <- flowDensity(x, c("B750-A", "SSC-A"), position = c(T, NA), gates = c(tp_ct_19, NA))
  tp_cp_19_nonB <- flowDensity(x, c("B750-A", "SSC-A"), position = c(F, NA), gates = c(tp_ct_19, NA))


  l.BnonB_pop <- list(
    tp_cp_19_B,
    tp_cp_19_nonB
  )

  names(l.BnonB_pop) <- c("B", "NonB")

  l.BnonB_gate <- f.generate_gate(l.BnonB_pop, gate.type = gate.type)

  l.BnonB_gate_cood <- f.extractGate_boundary(l.BnonB_gate)





  if (show.plots == T & nrow(x) >10) {
    par(mfrow = c(1, 4))
    flowPlot(x, channels = c("B750-A", "SSC-A"), title = identifier(x))
    lapply(l.BnonB_gate_cood[2], function(x) {
      lines(x, col = "blue", lty = 1, lwd = 2)
    })
    # abline(v = tp_ct_19, col = "blue")
  }

  return(l.BnonB_gate)
}











#' select B cells and NonB cell
#' NonB cells then can be subject to gdT cells then  CD2 and CD3 gating
#'  x<- gs_pop_get_data(gs, "Lymph")[[13]]
#'
f.NonT_gating <- function(x, gate.type = "r",
                          show.plots = T) {

  tp_ct_3 <- f.roughGate(x, "B810-A", tinypeak.removal = 1 / 10, alpha = 1 / 50)



  if (tp_ct_3 > 3) { # when there are a lot of B cells
    tp_ct_3 <- 2.5
  }


  tp_cp_3_T <- flowDensity(x, c("B810-A", "SSC-A"), position = c(T, NA), gates = c(tp_ct_3, NA))
  tp_cp_3_nonT <- flowDensity(x, c("B810-A", "SSC-A"), position = c(F, NA), gates = c(tp_ct_3, NA))


  l.TnonT_pop <- list(
    tp_cp_3_T,
    tp_cp_3_nonT
  )

  names(l.TnonT_pop) <- c("T", "NonT")

  l.TnonT_gate <- f.generate_gate(l.TnonT_pop, gate.type = gate.type)

  l.TnonT_gate_cood <- f.extractGate_boundary(l.TnonT_gate)





  if (show.plots == T) {
    par(mfrow = c(1, 4))
    flowPlot(x, channels = c("B810-A", "SSC-A"), title = identifier(x))
    lapply(l.TnonT_gate_cood[2], function(x) {
      lines(x, col = "blue", lty = 1, lwd = 2)
    })
    # abline(v = tp_ct_3, col = "blue")
  }

  return(l.TnonT_gate)
}



f.NonT_gating2 <- function(x, gate.type = "r",
                          show.plots = T) {
  
  sample_name <- identifier(x)
  print(sample_name)
  
  if (show.plots == T & nrow(x) >10) {
    par(mfrow = c(1, 4))
    flowPlot(x, channels = c("B810-A", "B750-A"), title = identifier(x), 
             xlim = c(0,4.5),
             ylim = c(0,4.5))
    # abline(v = tp_ct_3, col = "blue")
  }  
  
  cat("STEP 1: separate T and B on CD3\n")
  tp_ct_3 <- f.roughGate(x, "B810-A", tinypeak.removal = 1 / 10, alpha = 1 / 100)
  
  tp_cp_3_T <- flowDensity(x, c("B810-A", "B750-A"), position = c(T, NA), gates = c(tp_ct_3, NA))

  pks <- getPeaks(tp_cp_3_T, "B810-A",tinypeak.removal = 1 / 5)
  
  if(length(pks$Peaks) >1 ) {
    tp_ct_3 <-f.roughGate(tp_cp_3_T, "B810-A", tinypeak.removal = 1 / 10, alpha = 1 / 100, thrd = 2)
    
  }
  

  if (show.plots == T & nrow(x) >10) {
    flowPlot(x, channels = c("B810-A", "B750-A"), title = "TB cut")
    abline(v = tp_ct_3, col = "blue")
  }  
  
  cat("STEP 2: Set Gate on CD19 to capture CD19+ cells spreading into CD3 pos region\n")
  tp_cp_3_T <- flowDensity(x, c("B810-A", "B750-A"), position = c(T, NA), gates = c(tp_ct_3, NA))
  tp_ct_19 <- f.roughGate(tp_cp_3_T, "B750-A", tinypeak.removal = 1 / 10, alpha = 1 / 100, thrd = 1.5)
  
  if(is.infinite(tp_ct_19)){
    tp_ct_19 <- f.roughGate(x, "B750-A", tinypeak.removal = 1 / 10, alpha = 1 / 100, thrd = 1.5)
  }
  
  if (show.plots == T & nrow(x) >10) {
    flowPlot(x, channels = c("B810-A", "B750-A"), title = identifier(x))
    abline(v = tp_ct_3,h = tp_ct_19, col = "blue")
  }  
  
  print(tp_ct_3)
  cp_T <- flowDensity(x, c("B810-A", "B750-A"), position = c(T, F), gates = c(tp_ct_3, tp_ct_19))
  
  
  if(grepl("QC", sample_name, ignore.case = T) & grepl("CHEX", sample_name, ignore.case = T) ) {
    warning("QC Sample")
    cp_B <- flowDensity(x, c("B810-A", "B750-A"), position = c(F, T), gates = c(tp_ct_3, tp_ct_19))
    cp_nonTB <- flowDensity(x, c("B810-A", "B750-A"), position = c(F, F), gates = c(tp_ct_3, tp_ct_19))
    
    l.pop <- list(
      cp_B,
      cp_T,
      cp_nonTB
    )
    names(l.pop) <- c("B", "T", "NonTB")

  } else {
  
  l.pop <- list(
    cp_T,
    cp_T
  )
  names(l.pop) <- c("NonT", "T")
  
  }
  

  # l.pop_gate <- f.generate_gate(l.pop, gate.type = "r")
  
  
  l.pop_gate <- lapply(1:length(l.pop), function(i) { 
    g <- apply(l.pop[[i]]@filter, 2, range)
    #g[1, 2] <- -Inf
    #g[2,1] <- Inf
    rectangleGate(
      filterId = names(l.pop)[i],
      .gate = g
    )
  })
  names(l.pop_gate) <- names(l.pop)
  
  
  l.pop_gate_cood <- f.extractGate_boundary(l.pop_gate)
  
  
  if (show.plots == T & nrow(x) >10) {
    flowPlot(x, channels = c("B810-A", "B750-A"), title = identifier(x), xlim = c(-0.2, 4.5), ylim = c(-0.2, 4.5))
    lapply(l.pop_gate_cood, function(x) {
      lines(x, col = "blue", lty = 1, lwd = 2)
    })
  }  
  
  
  return(l.pop_gate)
  
  
}







#x <- cytoframe_to_flowFrame( fs[[1]]); gate.type = "q"; show.plots = T

f.Bgate <- function(x, gate.type = "r",
                    show.plots = T) {
  # plotDens(x, c("B750-A", "SSC-A"), density.overlay = c(T,T)  )
  
  
  
  ct_CD3 <-f.roughGate(x, "B810-A", thrd = 2, tinypeak.removal = 1 / 10, alpha = 1 / 50) #CD3
  ct_CD19 <- f.roughGate(x, "B750-A", thrd = 2, tinypeak.removal = 1 / 10, alpha = 1 / 50) # CD19
  
  
  tp_cp_CD19_NEG_CD3_POS <- flowDensity(x, c("B750-A", "B810-A"), position = c(F, T), gates = c(ct_CD19, ct_CD3))
  
  
  
  
  ct_CD20 <- f.roughGate(tp_cp_CD19_NEG_CD3_POS, "UV736-A", thrd = 2, tinypeak.removal = 1 / 10, alpha = 1 / 50) # CD19
  ct_CD22 <- f.roughGate(tp_cp_CD19_NEG_CD3_POS, "V785-A", thrd = 1.8, tinypeak.removal = 1 / 10, alpha = 1 / 50) # CD19
  
  
  
  
  if (show.plots == T & nrow(x) >10 & tp_cp_CD19_NEG_CD3_POS@cell.count >10) {
    par(mfrow = c(1, 4))
    flowPlot(x, channels = c("B750-A", "B810-A"), title = identifier(x), xlim = c(0,4.5), ylim = c(0,4.5))
    abline(v = ct_CD19, h = ct_CD3, col = "blue")

    flowPlot(x, channels = c("B750-A", "UV736-A"), title = identifier(x), xlim = c(0,4.5), ylim = c(0,4.5))
    abline(v = ct_CD19, h = ct_CD20, col = "blue")
    flowPlot(x, channels = c("B750-A", "V785-A"), title = identifier(x), xlim = c(0,4.5), ylim = c(0,4.5))
    abline(v = ct_CD19, h = ct_CD22, col = "blue")
    flowPlot(tp_cp_CD19_NEG_CD3_POS, channels = c("B750-A", "V785-A"), title = identifier(x), xlim = c(0,4.5), ylim = c(0,4.5))
  }

  
  CD19p <- flowDensity(x, c("B750-A", "B810-A"), position = c(T, F), gates = c(ct_CD19, ct_CD3))
  CD20p <- flowDensity(x, c("B750-A", "UV736-A"), position = c(T, T), gates = c(ct_CD19, ct_CD20))
  CD22p <- flowDensity(x, c("B750-A", "V785-A"), position = c(T, T), gates = c(ct_CD19, ct_CD22))
  
  
  l.pop <- list(
    CD19p,
    CD20p,
    CD22p
  )
  names(l.pop) <- c("CD19+", "CD20+", "CD22+")
  
  
  l.pop_gate <- f.generate_gate(l.pop, gate.type = gate.type)
  
  
  l.pop_gate_cood <- f.extractGate_boundary(l.pop_gate)
  
  
  if (show.plots == T & nrow(x) >10 & tp_cp_CD19_NEG_CD3_POS@cell.count >10) {
    par(mfrow = c(1, 4))
    flowPlot(x, channels = c("B750-A", "B810-A"), title = identifier(x), xlim = c(0,4.5), ylim = c(0,4.5))
    lines(l.pop_gate_cood[[1]], col = "blue", lty = 1, lwd = 2)
    
    flowPlot(x, channels = c("B750-A", "UV736-A"), title = identifier(x), xlim = c(0,4.5), ylim = c(0,4.5))
    lines(l.pop_gate_cood[[2]], col = "blue", lty = 1, lwd = 2)    
    flowPlot(x, channels = c("B750-A", "V785-A"), title = identifier(x), xlim = c(0,4.5), ylim = c(0,4.5))
    lines(l.pop_gate_cood[[3]], col = "blue", lty = 1, lwd = 2)    
    flowPlot(tp_cp_CD19_NEG_CD3_POS, channels = c("B750-A", "V785-A"), title = identifier(x), xlim = c(0,4.5), ylim = c(0,4.5))
    #lines(l.pop_gate_cood[[4]], col = "blue", lty = 1, lwd = 2)
  }
  
  
  return(l.pop_gate)
  
  
}


#x <- cytoframe_to_flowFrame( fs[[1]]); gate.type = "q"; show.plots = T
f.CD3_CD19_gating <- function(x, gate.type = "r",
                              show.plots = T) {
  
  ct_CD3 <-f.roughGate(x, "B810-A", thrd = 2, tinypeak.removal = 1 / 10, alpha = 1 / 50) #CD3
  ct_CD19 <- f.roughGate(x, "B750-A", thrd = 2, tinypeak.removal = 1 / 10, alpha = 1 / 50) # CD19
  
  
  tp_cp_CD19_POS_CD3_NEG <- flowDensity(x, c("B810-A", "B750-A"), position = c(F, T), gates = c(ct_CD3, ct_CD19))
  
  tp_cp_CD19_NEG_CD3_POS <- flowDensity(x, c("B810-A", "B750-A"), position = c(T, F), gates = c(ct_CD3, ct_CD19))
  
  
  if (show.plots == T & nrow(x) >10) {
    par(mfrow = c(1, 4))
    flowPlot(x, channels = c("B810-A", "B750-A"), title = identifier(x), xlim = c(0,4.5), ylim = c(0,4.5))
    abline(h = ct_CD19, v = ct_CD3, col = "blue")
    
  }
  
  l.pop <- list(
    tp_cp_CD19_POS_CD3_NEG,
    tp_cp_CD19_NEG_CD3_POS
  )
  
  names(l.pop) <- c("CD19pCD3n", "CD19nCD3p")
  
  
  l.pop_gate <- f.generate_gate(l.pop, gate.type = gate.type)
  
  
  l.pop_gate_cood <- f.extractGate_boundary(l.pop_gate)
  
  
  if (show.plots == T & nrow(x) >10) {
    flowPlot(x, channels = c("B810-A", "B750-A"), title = identifier(x), xlim = c(0,4.5), ylim = c(0,4.5))
    lines(l.pop_gate_cood[[1]], col = "blue", lty = 1, lwd = 2)
    lines(l.pop_gate_cood[[2]], col = "red", lty = 1, lwd = 2)  

  }
  
  
  
  
  return(l.pop_gate)

  
}






#' CD3 CD2 gating
#' NonB cells then can be subject to gdT cells then  CD2 and CD3 gating
#x<- cytoframe_to_flowFrame(fs[[1]]); gate.type = "r"; show.plots =T
#'
f.CD2CD3_gating <- function(x, gate.type = "r",
                            show.plots = T) {
  
  sample_name <- identifier(x)
  print(sample_name)
  
  # if(grepl("QC", sample_name) & grepl("CHEX", sample_name) & grepl("103", sample_name)) {
  #   warning("No PDC/Baso gate for CD-Chex CD103 Plus QC Sample")
  #   
  #   return(NULL)
  # } else {
  
  if (show.plots == T & nrow(x) > 10) {
    par(mfrow = c(1, 4))
    flowPlot(x, channels = c("B810-A", "B710-A"), title = identifier(x), xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5))
  }
  
  
  # pks <- getPeaks(x, "B810-A", tinypeak.removal = 1 / 10)
  # if (sum(pks$Peaks < 2) == 0) {
  #   tp_ct_CD3 <- deGate(x, "B810-A", tinypeak.removal = 1 / 10, use.upper = T, upper = F, alpha = 1 / 15)
  # } else {
  #   tp_ct_CD3 <- deGate(x, "B810-A", tinypeak.removal = 1 / 10)
  # }
  # 
  
  
  tp_ct_CD3 <- f.roughGate(x, "B810-A", thrd = 2, tinypeak.removal = 1 / 10, alpha = 1 / 50)
  
  

  tp_cp_CD3_POS <- flowDensity(x, c("B810-A", "B710-A"), position = c(T, NA), gates = c(tp_ct_CD3, NA))
  tp_cp_CD3_NEG <- flowDensity(x, c("B810-A", "B710-A"), position = c(F, NA), gates = c(tp_ct_CD3, NA))
  # flowPlot(tp_cp_CD3_POS, c("B810-A", "B710-A"), density.overlay = T  )

  
  if (show.plots == T  & tp_cp_CD3_NEG@cell.count > 10 ) {

    flowPlot(tp_cp_CD3_NEG, channels = c("B810-A", "B710-A"))

  }

  
  
  # if (sum(getPeaks(tp_cp_CD3_POS, "B710-A")$Peaks < 2.5) == 0) { # if there is no CD3+CD2- pop 
  #   tp_ct_CD3_POS_CD2 <- deGate(tp_cp_CD3_POS, "B710-A", use.upper = T, upper = F) # first peak
  #   tp_ct_CD3_POS_CD2_neg <- deGate(tp_cp_CD3_POS, "B710-A", use.upper = T, upper = F, alpha = 1 / 200) # push the cut a bit further down
  # } else {
  #   tp_ct_CD3_POS_CD2 <- tp_ct_CD3_POS_CD2_neg <- deGate(tp_cp_CD3_POS, "B710-A")
  # }
  # 
  
  
  
  tp_ct_CD3_POS_CD2 <- f.roughGate(x, "B710-A", thrd = 2.5)
  
  
  tp_ct_CD3_POS_CD2_neg <- f.roughGate(x, "B710-A", thrd = 2.5, alpha = 1 / 200) # push the cut a bit further down
  


  tp_cp_CD3_POS_CD2_POS <- flowDensity(tp_cp_CD3_POS, c("B810-A", "B710-A"), position = c(NA, T), gates = c(NA, tp_ct_CD3_POS_CD2))
  tp_cp_CD3_POS_CD2_NEG <- flowDensity(tp_cp_CD3_POS, c("B810-A", "B710-A"), position = c(NA, F), gates = c(NA, tp_ct_CD3_POS_CD2_neg))
  # flowPlot(tp_cp_CD3_POS_CD2_POS, c("B810-A", "B710-A"))



  # pks.cd2 <- getPeaks(tp_cp_CD3_NEG, "B710-A", tinypeak.removal = 1 / 10)
  # if (sum(pks.cd2$Peaks < 3) == 0) {
  #   tp_ct_CD3_NEG_CD2 <- deGate(tp_cp_CD3_NEG, "B710-A", tinypeak.removal = 1 / 10, use.upper = T, upper = F)
  # } else {
  #   tp_ct_CD3_NEG_CD2 <- deGate(tp_cp_CD3_NEG, "B710-A", all.cuts = T)
  # }

  tp_ct_CD3_NEG_CD2 <- f.roughGate(tp_cp_CD3_NEG, "B710-A", thrd = 2.5)

  # tp_ct_CD3_NEG_CD2.rot <-  tp_ct_CD3_NEG_CD2.rot[which.min(abs(tp_ct_CD3_NEG_CD2.rot -3))]

  temp <- flowDensity(tp_cp_CD3_NEG, c("B810-A", "B710-A"), position = c(NA, T), gates = c(NA, tp_ct_CD3_NEG_CD2))
  # flowPlot(tp_cp_CD4_POS_CD8_Neg, c("B810-A", "YG730-A"))

  if(temp@cell.count > 20) {

      tp_ct_CD3_NEG_CD2_pos <- f.roughGate(temp, "B710-A", thrd = 2)
    
    
  }else{tp_ct_CD3_NEG_CD2_pos <- -Inf}

  tp_ct_CD3_NEG_CD2 <- f.nearest_to_reference(c(tp_ct_CD3_NEG_CD2, tp_ct_CD3_NEG_CD2_pos), 3)



  tp_cp_CD3_NEG_CD2_NEG <- flowDensity(tp_cp_CD3_NEG, c("B810-A", "B710-A"), position = c(NA, F), gates = c(NA, tp_ct_CD3_NEG_CD2))
  tp_cp_CD3_NEG_CD2_POS <- flowDensity(tp_cp_CD3_NEG, c("B810-A", "B710-A"), position = c(NA, T), gates = c(NA, tp_ct_CD3_NEG_CD2))

  # ROTATE BACK
  # tp_cp_CD3_NEG_CD2_NEG <- rotate.fd(tp_ct_CD3_NEG_CD2_NEG.rot, angle = rotate.theta)
  # tp_cp_CD3_NEG_CD2_POS <- rotate.fd(tp_ct_CD3_NEG_CD2_POS.rot, angle = rotate.theta)



  # flowPlot(tp_cp_CD4_POS_CD8_Neg, c("B537-A", "YG730-A"))
  # use inflection point to tight the CD4 gating

  l.2and3_pop <- list(
    tp_cp_CD3_POS_CD2_POS,
    tp_cp_CD3_POS_CD2_NEG,
    tp_cp_CD3_NEG_CD2_POS,
    tp_cp_CD3_NEG_CD2_NEG
  )

  names(l.2and3_pop) <- c("CD3pCD2p", "CD3pCD2n", "CD3nCD2p", "CD3nCD2n")

  # remove gates that has less than 3 cells

  cn <- unlist(lapply(l.2and3_pop, function(x) {
    x@cell.count
  }))

  # l.2and3_pop <- l.2and3_pop[cn>3]


  l.2and3_gate <- f.generate_gate(l.2and3_pop, gate.type = gate.type)
  l.2and3_gate <- l.2and3_gate
  l.2and3_gate_cood <- f.extractGate_boundary(l.2and3_gate[cn > 3])



  if (show.plots == T & nrow(x) >10) {
    flowPlot(x, channels = c("B810-A", "B710-A"), xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5))
    lapply(l.2and3_gate_cood, function(x) {
      lines(x, col = "blue", lty = 1, lwd = 2)
    })
  }

  return(l.2and3_gate)
}



#'  x<- gs_pop_get_data(gs, "abT")[[1]]
#'
f.CD4CD8_gating <- function(x, gate.type = "r",
                            show.plots = T) {
  # 4 and 8 is comparatively easy to gate
  # pks.x <- getPeaks(x, "B537-A")
  # pks.y <- getPeaks(x, "YG730-A")


  # if (sum(pks.y$Peaks > 0) == 1) {
  #   tp_ct_CD4 <- deGate(x, "YG730-A", use.upper = T, upper = F)
  # } else {
  #   tp_ct_CD4 <- deGate(x, "YG730-A", all.cuts = T)
  # }
  # tp_ct_CD4 <- f.nearest_to_reference(tp_ct_CD4, 2.5)
  # 
  
  tp_ct_CD4 <- f.roughGate(x, "YG730-A", thrd = 2)

  
  tp_ct_CD8 <- f.roughGate(x, "B537-A", thrd = 2)



  # if (sum(pks.x$Peaks > 0) == 1) {
  #   tp_ct_CD8 <- deGate(x, "B537-A", use.upper = T, upper = F)
  # } else {
  #   tp_ct_CD8 <- deGate(x, "B537-A", all.cuts = T)
  # }
  # 
  # 
  # 
  # tp_ct_CD8 <- f.nearest_to_reference(tp_ct_CD8, 2)
  # 


  l.4and8_pop <- lapply(list(c(T, T), c(F, T), c(T, F), c(F, F)), function(p) {
    flowDensity(x, c("B537-A", "YG730-A"), position = p, gates = c(tp_ct_CD8, tp_ct_CD4))
  })
  names(l.4and8_pop) <- c("CD4pCD8p", "CD4pCD8n", "CD4nCD8p", "CD4nCD8n")


  # fine tuning CD4 Tpop
  tp_cp_CD4_POS_CD8_Neg <- l.4and8_pop[[2]]
  tp_ct_D4_POS_CD8_Neg.CD4 <- f.roughGate(tp_cp_CD4_POS_CD8_Neg, thrd = 2, "YG730-A", alpha = 1 / 50)


  tp_cp_CD4_Neg_CD8_Pos <- l.4and8_pop[[3]]
  tp_ct_D4_NEG_CD8Pos.CD4 <- f.roughGate(tp_cp_CD4_Neg_CD8_Pos,thrd = 2, "YG730-A", alpha = 1 / 30)
  
  
  if(!is.infinite(tp_ct_D4_POS_CD8_Neg.CD4)) { 
    l.4and8_pop[[2]] <- flowDensity(tp_cp_CD4_POS_CD8_Neg, c("B537-A", "YG730-A"),
                                    position = c(NA, T), gates = c(NA, tp_ct_D4_POS_CD8_Neg.CD4)
    )
  }

  if(!is.infinite(tp_ct_D4_NEG_CD8Pos.CD4)) {
    l.4and8_pop[[3]] <- flowDensity(tp_cp_CD4_Neg_CD8_Pos, c("B537-A", "YG730-A"),
                                    position = c(NA, F), gates = c(NA, tp_ct_D4_NEG_CD8Pos.CD4)
    )
  }
  

  l.4and8_gate <- f.generate_gate(l.4and8_pop, gate.type = gate.type)
  l.4and8_gate_cood <- f.extractGate_boundary(l.4and8_gate)

  if (show.plots == T & nrow(x) > 2) {
    par(mfrow = c(1, 4))
    flowPlot(x, channels = c("B537-A", "YG730-A"), xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5), title = identifier(x))
    abline(v = tp_ct_CD8, h = tp_ct_CD4, col = "green")
    lapply(l.4and8_gate_cood, function(x) {
      lines(x, col = "blue", lty = 1, lwd = 2)
    })
  }

  return(l.4and8_gate)
}




f.CD4_LHES_gating <- function(x, gate.type = "r", gate.name, x.axis = x.axis, default_ct_CD4_high = 2.5, 
                            show.plots = T) {

  tp_ct_CD4 <- f.roughGate(x, "YG730-A", thrd = 2, alpha = 1/100)
  tp_ct_X <- f.roughGate(x, x.axis, thrd = 2)
  
  tp_cp_CD4_POS_X_Neg <- flowDensity(x, c(x.axis, "YG730-A"), position = c(NA,T), gates = c(tp_ct_X, tp_ct_CD4))

  
  
  
  if (show.plots == T & nrow(x) > 10) {
    par(mfrow = c(1, 4))
    flowPlot(x, channels = c(x.axis, "YG730-A"), xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5), title = identifier(x))
     abline(v = tp_ct_X, h = tp_ct_CD4, col = "green")

  }
  
  tp_cp_CD4_POS_X_POS <- flowDensity(x, c(x.axis, "YG730-A"), position = c(T,T), gates = c(tp_ct_X, tp_ct_CD4))
  #tp_ct_CD4_thrd <- f.roughGate(tp_cp_CD4_POS_X_POS, "YG730-A", thrd = default_ct_CD4_high)
  #cat("CD4 threshold from X+CD4+ pop:", tp_ct_CD4_thrd, "\n")
  #tp_ct_CD4_high <- f.roughGate(tp_cp_CD4_POS_X_Neg, "YG730-A", thrd = max(tp_ct_CD4_thrd, default_ct_CD4_high))
  
  
  if (show.plots == T & tp_cp_CD4_POS_X_Neg@cell.count > 10) {

    flowPlot(tp_cp_CD4_POS_X_Neg, channels = c(x.axis, "YG730-A"), xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5), title = identifier(x))
    abline(v = tp_ct_X, h = default_ct_CD4_high, col = "green")
  } else {
    plot(NULL, xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5), xlab = "", ylab = "", main = "Not enough cells")
  }
  
  cp_CD4highpop <- flowDensity(x, c(x.axis, "YG730-A"), position = c(F,T), gates = c(tp_ct_X, max(default_ct_CD4_high, tp_ct_CD4)))
  
  
  l.pop = list(cp_CD4highpop)
  names(l.pop) <- gate.name
  
  
  l.gate <- f.generate_gate(l.pop, gate.type = gate.type)
  l.gate_cood <- f.extractGate_boundary(l.gate)
  
  
  if (show.plots == T & nrow(x) > 10) {
    flowPlot(x, channels = c(x.axis, "YG730-A"), xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5), title = identifier(x))
    lapply(l.gate_cood, function(x) {
      lines(x, col = "blue", lty = 1, lwd = 2)
    })
  }
  
  return(l.gate)
 # return(l.4and8_gate)
}




f.CD4CD7_LHES_gating <- function(x, gate.type = "r",gate.name,
                                 show.plots = T) {
  
  tp_ct_CD4 <- f.roughGate(x, "YG730-A", thrd = 2)
  tp_ct_CD7 <- f.roughGate(x, "UV736-A", thrd = 2)
  
  tp_cp_CD4_POS_CD7_Neg <- flowDensity(x, c("UV736-A", "YG730-A"), position = c(F,T), gates = c(tp_ct_CD7, tp_ct_CD4))
  
  
  
  if (show.plots == T & nrow(x) > 10) {
    par(mfrow = c(1, 4))
    flowPlot(x, channels = c("UV736-A", "YG730-A"), xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5), title = identifier(x))
    abline(v = tp_ct_CD7, h = tp_ct_CD4, col = "green")
    
  }
  
  tp_ct_CD4_high <- f.roughGate(tp_cp_CD4_POS_CD7_Neg, "YG730-A", thrd = 3, tinypeak.removal = 1/10)
  
  
  if (show.plots == T & tp_cp_CD4_POS_CD7_Neg@cell.count > 10) {
    
    flowPlot(tp_cp_CD4_POS_CD7_Neg, channels = c("UV736-A", "YG730-A"), xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5), title = identifier(x))
    abline(v = tp_ct_CD7, h = tp_ct_CD4_high, col = "green")
  }
  
  cp_CD4highpop <- flowDensity(x, c("UV736-A", "YG730-A"), position = c(F,T), gates = c(tp_ct_CD7,  max(tp_ct_CD4_high, tp_ct_CD4)))
  
  
  l.pop = list(cp_CD4highpop)
  names(l.pop) <- gate.name
  
  l.gate <- f.generate_gate(l.pop, gate.type = gate.type)
  l.gate_cood <- f.extractGate_boundary(l.gate)
  
  
  if (show.plots == T & nrow(x) > 10) {
    flowPlot(x, channels = c("UV736-A", "YG730-A"), xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5), title = identifier(x))
    lapply(l.gate_cood, function(x) {
      lines(x, col = "blue", lty = 1, lwd = 2)
    })
  }
  
  return(l.gate)
  # return(l.4and8_gate)
}





f.CD4CD5_LHES_gating <- function(x, gate.type = "r",gate.name,
                                 show.plots = T) {
  
  tp_ct_CD4 <- f.roughGate(x, "YG730-A", thrd = 2)
  tp_ct_CD5 <- f.roughGate(x, "UV809-A", thrd = 2)
  
  tp_cp_CD4_POS_CD5_Neg <- flowDensity(x, c("UV809-A", "YG730-A"), position = c(F,T), gates = c(tp_ct_CD5, tp_ct_CD4))
  
  
  
  if (show.plots == T & nrow(x) > 10) {
    par(mfrow = c(1, 4))
    flowPlot(x, channels = c("UV809-A", "YG730-A"), xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5), title = identifier(x))
    abline(v = tp_ct_CD5, h = tp_ct_CD4, col = "green")
    
  }
  
  tp_ct_CD4_high <- f.roughGate(tp_cp_CD4_POS_CD5_Neg, "YG730-A", thrd = 3, tinypeak.removal = 1/10)
  
  
  if (show.plots == T & tp_cp_CD4_POS_CD5_Neg@cell.count > 10) {
    
    flowPlot(tp_cp_CD4_POS_CD5_Neg, channels = c("UV809-A", "YG730-A"), xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5), title = identifier(x))
    abline(v = tp_ct_CD5, h = tp_ct_CD4_high, col = "green")
  }
  
  cp_CD4highpop <- flowDensity(x, c("UV809-A", "YG730-A"), position = c(F,T), gates = c(tp_ct_CD5,  max(tp_ct_CD4_high, tp_ct_CD4)))
  
  
  l.pop = list(cp_CD4highpop)
  names(l.pop) <- gate.name
  
  l.gate <- f.generate_gate(l.pop, gate.type = gate.type)
  l.gate_cood <- f.extractGate_boundary(l.gate)
  
  
  if (show.plots == T & nrow(x) > 10) {
    flowPlot(x, channels = c("UV809-A", "YG730-A"), xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5), title = identifier(x))
    lapply(l.gate_cood, function(x) {
      lines(x, col = "blue", lty = 1, lwd = 2)
    })
  }
  
  return(l.gate)
  # return(l.4and8_gate)
}




f.CD4CD14_LHES_gating <- function(x, gate.type = "r",gate.name, default_ct_CD4 = 2, 
                                 show.plots = T) {
  
  tp_ct_CD4 <- f.roughGate(x, "YG730-A", thrd = 2)
  tp_ct_CD14 <- f.roughGate(x, "B602-A", thrd = 2)
  
  tp_cp_CD4_POS_CD14_Neg <- flowDensity(x, c("B602-A", "YG730-A"), position = c(F,T), gates = c(tp_ct_CD14, tp_ct_CD4))
  
  
  
  if (show.plots == T & nrow(x) > 10) {
    par(mfrow = c(1, 4))
    flowPlot(x, channels = c("B602-A", "YG730-A"), xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5), title = identifier(x))
    abline(v = tp_ct_CD14, h = tp_ct_CD4, col = "green")
    
  }
  
  tp_ct_CD4_high <- f.roughGate(tp_cp_CD4_POS_CD14_Neg, "YG730-A", thrd = default_ct_CD4, tinypeak.removal = 1/10)
  
  
  if (show.plots == T & tp_cp_CD4_POS_CD14_Neg@cell.count > 10) {
    
    flowPlot(tp_cp_CD4_POS_CD14_Neg, channels = c("B602-A", "YG730-A"), xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5), title = identifier(x))
    abline(v = tp_ct_CD14, h = tp_ct_CD4_high, col = "green")
  }
  
  cp_CD4highpop <- flowDensity(x, c("B602-A", "YG730-A"), position = c(F,T), gates = c(tp_ct_CD14,  max(tp_ct_CD4_high, tp_ct_CD4)))
  
  
  l.pop = list(cp_CD4highpop)
  names(l.pop) <- gate.name
  
  l.gate <- f.generate_gate(l.pop, gate.type = gate.type)
  l.gate_cood <- f.extractGate_boundary(l.gate)
  
  
  if (show.plots == T & nrow(x) > 10) {
    flowPlot(x, channels = c("B602-A", "YG730-A"), xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5), title = identifier(x))
    lapply(l.gate_cood, function(x) {
      lines(x, col = "blue", lty = 1, lwd = 2)
    })
  }
  
  return(l.gate)
  # return(l.4and8_gate)
}






f.CD4CD26_CTCL_gating <- function(x, gate.type = "r",gate.name,
                                  show.plots = T) {
  
  tp_ct_CD4 <- f.roughGate(x, "YG730-A", thrd = 3)
  tp_ct_CD26 <- f.roughGate(x, "YG602-A", thrd = 2)
  
  tp_cp_CD4_POS_CD26_Neg <- flowDensity(x, c("YG730-A", "YG602-A"), position = c(T,F), gates = c(tp_ct_CD4, tp_ct_CD26))
  
  
  
  if (show.plots == T & nrow(x) > 10) {
    par(mfrow = c(1, 4))
    flowPlot(x, channels = c("YG730-A", "YG602-A"), xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5), title = identifier(x))
    abline(v = tp_ct_CD4, h = tp_ct_CD26, col = "green")
    
  }
  
  #tp_ct_CD4_high <- f.roughGate(tp_cp_CD4_POS_CD14_Neg, "YG730-A", thrd = 3, tinypeak.removal = 1/10)
  
  
  # if (show.plots == T & tp_cp_CD4_POS_CD14_Neg@cell.count > 10) {
  #   
  #   flowPlot(tp_cp_CD4_POS_CD14_Neg, channels = c("YG730-A", "YG602-A"), xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5), title = identifier(x))
  #   abline(v = tp_ct_CD26, h = tp_ct_CD4_high, col = "green")
  # }
  # 
  # cp_CD4highpop <- flowDensity(x, c("YG730-A", "YG602-A"), position = c(F,T), gates = c(tp_ct_CD26, tp_ct_CD4_high))
  # 
  
  l.pop = list(tp_cp_CD4_POS_CD26_Neg)
  names(l.pop) <- gate.name
  
  l.gate <- f.generate_gate(l.pop, gate.type = gate.type)
  l.gate_cood <- f.extractGate_boundary(l.gate)
  
  
  if (show.plots == T & nrow(x) > 10) {
    flowPlot(x, channels = c("YG730-A", "YG602-A"), xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5), title = identifier(x))
    lapply(l.gate_cood, function(x) {
      lines(x, col = "blue", lty = 1, lwd = 2)
    })
  }
  
  return(l.gate)
  # return(l.4and8_gate)
}








#' NK CD56 gating
#' NonB cells then can be subject to gdT cells then  CD2 and X gating
#'  x<- gs_pop_get_data(gs, "NonB")[[1]]
#' default_ct_CD3 <- l.ct_CD3[[1]]
f.CD3_CD56_gating <- function(x, gate.type = "r", default_ct_CD3,
                              show.plots = T) {
  print(identifier(x))
  if (show.plots == T & nrow(x) > 10) {
    par(mfrow = c(1, 4))
    flowPlot(x, channels = c("B810-A", "V470-A"), xlim = c(0, 4.5), ylim = c(0, 4.5))

  }
  

  pks <- getPeaks(x, "V470-A", tinypeak.removal = 1 / 100)
  if (length(pks$Peaks) == 1) { # inflection gate
    tp_ct_CD56 <- deGate(x, "V470-A", tinypeak.removal = 1 / 100, use.upper = T, upper = T)
  } else {
    tp_ct_CD56 <- deGate(x, "V470-A", tinypeak.removal = 1 / 100, all.cuts = T)
  }

  if (missing(default_ct_CD3)) {
    tp_ct_CD3 <- deGate(x, "B810-A")
  } else {
    tp_ct_CD3 <- default_ct_CD3
  }

  tp_cp_CD3_NEG <- flowDensity(x, c("B810-A", "V470-A"), position = c(F, NA), gates = c(tp_ct_CD3, NA))
  # plotDens(tp_cp_CD3_NEG, c("B810-A", "V470-A"), density.overlay = c(T,T)  )



  # rotate.theta <- 0 # just to tilt a bit to avoid the diagnoal pops
  # tp_cp_CD3CD56_NEG.rot <- rotate.data(getflowFrame(tp_cp_CD3_NEG), c("B810-A", "V470-A"), theta = rotate.theta)$data
  # 
  # pks <- getPeaks(tp_cp_CD3_NEG, "V470-A", tinypeak.removal = 1 / 100)
  # if (length(pks$Peaks) == 1) { # inflection gate
  #   tp_ct_CD3_NEG_CD56 <- deGate(tp_cp_CD3_NEG, "V470-A", tinypeak.removal = 1 / 100, use.upper = T, upper = T)
  # } else {
  #   tp_ct_CD3_NEG_CD56 <- deGate(tp_cp_CD3_NEG, "V470-A", tinypeak.removal = 1 / 100, all.cuts = T)
  # }
  # 
  # ct_CD56_V470_NK_preSet <- 1
  # 
  # tp_ct_CD3_NEG_CD56 <- f.nearest_to_reference(c(tp_ct_CD56, tp_ct_CD3_NEG_CD56), ct_CD56_V470_NK_preSet)

  
  
  tp_ct_CD3_NEG_CD56 <- f.roughGate(tp_cp_CD3_NEG, "V470-A", thrd = 1.5, alpha = 1/2)
  # tp_ct_CD3_NEG_CD56[which.min(abs(tp_ct_CD3_NEG_CD56 - ct_CD56_V470_NK_preSet))]

  tp_cp_CD3_NEG_CD56_POS <- flowDensity(tp_cp_CD3_NEG, c("B810-A", "V470-A"), position = c(NA, T), gates = c(NA, tp_ct_CD3_NEG_CD56))
  tp_cp_CD3_POS_CD56_POS <- flowDensity(x, c("B810-A", "V470-A"), position = c(T, T), gates = c(tp_ct_CD3, tp_ct_CD3_NEG_CD56))
  

  l.56and3_gate <- f.generate_gate(list(`CD56pCD3n` = tp_cp_CD3_NEG_CD56_POS,
                                        `NKT` = tp_cp_CD3_POS_CD56_POS
                                        ), gate.type = gate.type)
  # l.56and3_gate$NK@max[2] <- 4.3
  l.56and3_gate_cood <- f.extractGate_boundary(l.56and3_gate)



  if (show.plots == T & nrow(x) > 10 & tp_cp_CD3_NEG@cell.count > 10 ) {
    flowPlot(tp_cp_CD3_NEG, channels = c("B810-A", "V470-A"), xlim = c(0, 4.5), ylim = c(0, 4.5))
    flowPlot(x, channels = c("B810-A", "V470-A"), xlim = c(0, 4.5), ylim = c(0, 4.5))
    lapply(l.56and3_gate_cood, function(x) {
      lines(x, col = "blue", lty = 1, lwd = 2)
    })
  }

  return(l.56and3_gate)
}


#' NK CD56 gating
#' NonB cells then can be subject to gdT cells then  CD2 and CD3 gating
#'  x<- gs_pop_get_data(gs, "NonB")[[1]]
#' default_ct_CD3 <- l.ct_CD3[[1]]
f.CD3_CD16_gating <- function(x, gate.type = "r", default_ct_CD3,
                              show.plots = T) {
  print(identifier(x))
  # flowPlot(x, c("B810-A", "V470-A"), density.overlay = c(T,T)  )
  
  pks <- getPeaks(x, "V710-A", tinypeak.removal = 1 / 100)
  if (length(pks$Peaks) == 1) { # inflection gate
    tp_ct_CD16 <- deGate(x, "V710-A", tinypeak.removal = 1 / 100, use.upper = T, upper = T)
  } else {
    tp_ct_CD16 <- deGate(x, "V710-A", tinypeak.removal = 1 / 100, all.cuts = T)
  }
  
  if (missing(default_ct_CD3)) {
    tp_ct_CD3 <- deGate(x, "B810-A")
  } else {
    tp_ct_CD3 <- default_ct_CD3
  }
  
  tp_cp_CD3_NEG <- flowDensity(x, c("B810-A", "V710-A"), position = c(F, NA), gates = c(tp_ct_CD3, NA))
  # plotDens(tp_cp_CD3_NEG, c("B810-A", "V710-A"), density.overlay = c(T,T)  )
  
  
  
  # plotDens(tp_cp_CD3_NEG, c("B810-A", "B710-A"), density.overlay = c(T,T)  )
  rotate.theta <- 0 # just to tilt a bit to avoid the diagnoal pops
  tp_cp_CD3CD16_NEG.rot <- rotate.data(getflowFrame(tp_cp_CD3_NEG), c("B810-A", "V710-A"), theta = rotate.theta)$data
  # flowPlot(tp_cp_CD3CD16_NEG.rot, c("B810-A", "V710-A"), density.overlay = c(T,T)  )
  
  
  
  pks <- getPeaks(tp_cp_CD3_NEG, "V710-A", tinypeak.removal = 1 / 100)
  if (length(pks$Peaks) == 1) { # inflection gate
    tp_ct_CD3_NEG_CD16 <- deGate(tp_cp_CD3_NEG, "V710-A", tinypeak.removal = 1 / 100, use.upper = T, upper = T)
  } else {
    tp_ct_CD3_NEG_CD16 <- deGate(tp_cp_CD3_NEG, "V710-A", tinypeak.removal = 1 / 100, all.cuts = T)
  }
  
  ct_CD16_V710_NK_preSet <- 1.8
  
  tp_ct_CD3_NEG_CD16 <- f.nearest_to_reference(c(tp_ct_CD16, tp_ct_CD3_NEG_CD16), ct_CD16_V710_NK_preSet)
  
  # tp_ct_CD3_NEG_CD16[which.min(abs(tp_ct_CD3_NEG_CD16 - ct_CD16_V710_NK_preSet))]
  
  tp_cp_CD3_NEG_CD16_POS <- flowDensity(tp_cp_CD3_NEG, c("B810-A", "V710-A"), position = c(NA, T), gates = c(NA, tp_ct_CD3_NEG_CD16))
  
  
  l.16and3_gate <- f.generate_gate(list(`CD16pCD3n` = tp_cp_CD3_NEG_CD16_POS), gate.type = gate.type)
  # l.16and3_gate$NK@max[2] <- 4.3
  l.16and3_gate_cood <- f.extractGate_boundary(l.16and3_gate)
  
  
  
  if (show.plots == T & nrow(x) > 10 & tp_cp_CD3_NEG@cell.count > 10) {
    par(mfrow = c(1, 4))
    flowPlot(tp_cp_CD3_NEG, channels = c("B810-A", "V710-A"), xlim = c(-1, 4.5), ylim = c(0.5, 4.5))
    flowPlot(x, channels = c("B810-A", "V710-A"), xlim = c(-1, 4.5), ylim = c(0.5, 4.5))
    lapply(l.16and3_gate_cood, function(x) {
      lines(x, col = "blue", lty = 1, lwd = 2)
    })
  }
  
  return(l.16and3_gate)
}



# x<- gs_pop_get_data(gs, "CD3pCD2p")[[22]]
# default_ct_CD57 <- l.ct_CD57[[22]]
# axis.x <- "B510-A"
# axis.y <- "UV585-A"

f.quad_gating <- function(x, parent,
                          axis.x,
                          axis.y,
                          default_ct_x,
                          default_ct_y,
                          show.plots = T) {
  # flowPlot(x, c(axis.x, axis.y) )
  if (missing(default_ct_x)) {
    tp_ct_x <- deGate(x, axis.x)
  } else {
    tp_ct_x <- default_ct_x
  }

  if (missing(default_ct_y)) {
    # check where here is major bimodal ditribution
    pks <- getPeaks(x, axis.y, tinypeak.removal = 1 / 10)
    if (length(pks$Peaks) == 1) {
      tp_ct_y <- deGate(x, axis.y, tinypeak.removal = 1 / 10, use.upper = F, upper = T)
    } else {
      tp_ct_y <- deGate(x, axis.y)
    }
  } else {
    tp_ct_y <- default_ct_y
  }

  qg <- c(tp_ct_x, tp_ct_y)
  names(qg) <- c(axis.x, axis.y)

  l.gate <- list(quadGate(.gate = qg))
  names(l.gate) <- parent

  # l.56and3_gate$NK@max[2] <- 4.3
  # l.56and3_gate_cood <- f.extractGate_boundary(l.56and3_gate)


  if (show.plots == T) {
    par(mfrow = c(1, 4))
    flowPlot(x, channels = c(axis.x, axis.y), xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5))
    abline(v = tp_ct_x, h = tp_ct_y, col = "blue", lty = 1, lwd = 2)
  }

  return(l.gate)
}


# General gate
# debug_case <- 12
# x <- gs_pop_get_data(gs, "abT")[[debug_case]]
# default_ct_y <- l.ct_CD8[[debug_case]]
# gate.type ="r"
# axis.x <- "R675-A"
# axis.y <- "V750-A"
# mid.x <-1.5
# mid.y <- 2
# subpop <- list('CD8pCD57p' = c(T,T), `CD4pCD57p` = c(T,F))

f.general_gating <- function(x, gate.type = "r",
                             # parent,
                             children_gates,
                             axis.x, mid.x = 2.2, alpha.x = 1 / 10,
                             axis.y, mid.y = 2.2, alpha.y = 1 / 10,
                             tiny_peak = 1 / 25,
                             default_ct_x,
                             default_ct_y,
                             twin.factor =0.98,
                             show.plots = T) {
  # flowPlot(x, c(axis.x, axis.y) )
  initial_tiny_peak <- tiny_peak
  
  dat <- exprs(x)[,c(axis.x, axis.y)]
  
  exprs(x)[,c(axis.x, axis.y)] <- apply(dat, 2, function(x) {
    x <- pmax(x, 0) 
    x
  })
  

  twin.factor.x <- twin.factor
  twin.factor.y <- twin.factor
  
  if (show.plots == T & nrow(x) >10) {
    par(mfrow = c(1, 4))
    flowPlot(x,
             channels = c(axis.x, axis.y), # xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5),
             title = identifier(x)
    )
  }
    

  cat("Perform gate on x axis: ",axis.x,  "\n")
  
  if (missing(default_ct_x)) {
    tp_ct_x <- f.roughGate(x, axis.x, thrd = mid.x, tinypeak.removal = initial_tiny_peak, twin.factor = twin.factor.x, alpha = alpha.x)
  } else {
    cat("using default cut for x axis: ",axis.x,  "\n")
    tp_ct_x <- default_ct_x
  } 
  
  cat("Perform gate on y axis: ",axis.y,  "\n")
  
  if (missing(default_ct_y)) {
    tp_ct_y <-f.roughGate(x, axis.y, thrd = mid.y, tinypeak.removal = initial_tiny_peak, twin.factor = twin.factor.y, alpha = alpha.y)
  } else {
    cat("using default cut for y axis: ",axis.y,  "\n")
    tp_ct_y <- default_ct_y
  }
  

  
  if (missing(default_ct_x) &(!is.infinite(tp_ct_y)) & (tp_ct_x < 1/2*mid.x | tp_ct_x > 2*mid.x) )  {    
    cat("cutting point on x axis: ", axis.x, " is too far away\n")
    #tp_ct_x <- deGate(x, axis.x, n.sd = 0.5, sd.threshold = T)
    #tp_ct_x <- deGate(x, axis.x, n.sd = 1, sd.threshold = T)
    
    temp_y_pos <- flowDensity(x, c(axis.x, axis.y), position = c(NA, T), gates = c(NA, tp_ct_y))
    temp_y_neg <- flowDensity(x, c(axis.x, axis.y), position = c(NA, F), gates = c(NA, tp_ct_y))
    
    tp_ct_x_1 <- f.roughGate(temp_y_pos, axis.x, thrd = mid.x, tinypeak.removal = initial_tiny_peak, twin.factor = twin.factor.x, alpha = alpha.x)
    cat("x cuts from upper:", tp_ct_x_1, "\n")
    tp_ct_x_2 <- f.roughGate(temp_y_neg, axis.x, thrd = mid.x, tinypeak.removal = initial_tiny_peak, twin.factor = twin.factor.x, alpha = alpha.x)
    cat("x cuts from lower:", tp_ct_x_2, "\n")
    
    
    tp_ct_x_all <- c(tp_ct_x_1, tp_ct_x_2)
    tp_ct_x <- f.nearest_to_reference(tp_ct_x_all, mid.x)
    # flowPlot(temp_y, c(axis.x, axis.y) )
  }
  
  

  
  if(missing(default_ct_y) &(!is.infinite(tp_ct_y))& (tp_ct_y < 1/2*mid.y | tp_ct_y > 2*mid.y) ) {
    cat("cutting point on y axis: ", axis.y, " is too far away\n")    #tp_ct_y <- deGate(x, axis.y, n.sd = 0.5, sd.threshold = T)
    
    
    temp_x_pos <- flowDensity(x, c(axis.x, axis.y), position = c(T, NA), gates = c(NA, tp_ct_x))
    temp_x_neg <- flowDensity(x, c(axis.x, axis.y), position = c(F, NA), gates = c(NA, tp_ct_x))
    
    tp_ct_y_1 <- f.roughGate(temp_x_pos, axis.y, thrd = mid.y, tinypeak.removal = initial_tiny_peak, twin.factor = twin.factor.y, alpha = alpha.y)
    cat("x cuts from right: ", tp_ct_y_1, "\n")
    
    
    tp_ct_y_2 <- f.roughGate(temp_x_neg, axis.y, thrd = mid.y, tinypeak.removal = initial_tiny_peak, twin.factor = twin.factor.y, alpha = alpha.y)
    cat("x cuts from left:", tp_ct_y_2, "\n")
    tp_ct_y_all <- c(tp_ct_y_1, tp_ct_y_2)
    tp_ct_y <- f.nearest_to_reference(tp_ct_y_all, mid.y)
    
    
    # flowPlot(temp_y, c(axis.x, axis.y) )
  }
  
  
  
  # flowPlot(x, c(axis.x, axis.y) )
  # abline(h = tp_ct_y, v = tp_ct_x)


  if (gate.type == "q") { 
    qg <- c(tp_ct_x, tp_ct_y)
    names(qg) <- c(axis.x, axis.y)

    l.gate <- list(quadGate(.gate = qg)) 
    # names(l.gate) <- parent
  } else {
    l.pop <- lapply(children_gates, function(p) {
      flowDensity(x, c(axis.x, axis.y), position = p, gates = c(tp_ct_x, tp_ct_y))
    })
    names(l.pop) <- names(children_gates)
    l.gate <- f.generate_gate(l.pop, gate.type = gate.type)
    l.gate_cood <- f.extractGate_boundary(l.gate)
  }

  Gate_colors <- c("blue", "green", "magenta", "brown")
  if (show.plots == T & nrow(x) >10) {
    flowPlot(x,
      channels = c(axis.x, axis.y),  xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5),
      title = identifier(x)
    )
    abline(v = tp_ct_x, h = tp_ct_y, col = "black", lty = 1, lwd = 2)
    if (gate.type == "r" | gate.type == "p") {
      lapply(seq(length(l.gate_cood)), function(x) {
        lines(l.gate_cood[[x]],
          col = Gate_colors[[x]],
          lty = 1, lwd = 2
        )
      })
    }
  }

  return(l.gate)
}
















#x<- cytoframe_to_flowFrame(fs[[4]])

f.TRBC_gating <- function(x, gate.type = "r",
                             # parent,
                             children_gates,
                             axis.x, mid.x = 1.9, alpha.x = 1 / 10,
                             axis.y, mid.y = 2, alpha.y = 1 / 10,
                             tiny_peak = 1 / 25,
                             default_ct_x,
                             default_ct_y,
                             twin.factor =0.98,
                             show.plots = T) {
  # flowPlot(x, c(axis.x, axis.y) )
  initial_tiny_peak <- tiny_peak
  
  dat <- exprs(x)[,c(axis.x, axis.y)]
  
  exprs(x)[,c(axis.x, axis.y)] <- apply(dat, 2, function(x) {
    x <- pmax(x, 0) 
    x
  })
  
  
  twin.factor.x <- twin.factor
  twin.factor.y <- twin.factor
  
  if (show.plots == T & nrow(dat) > 10) {
    par(mfrow = c(1, 4))
    flowPlot(x,
             channels = c(axis.x, axis.y), # xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5),
             title = identifier(x)
    )
  }
  
  
  cat("Perform gate for ",axis.x,  "\n")
  
  if (missing(default_ct_x)) {
    tp_ct_x <- f.roughGate(x, axis.x, thrd = mid.x, tinypeak.removal = initial_tiny_peak, twin.factor = twin.factor.x, alpha = alpha.x)
  } else {
    cat("using default cut for",axis.x,  "\n")
    tp_ct_x <- default_ct_x
  } 
  
  cat("Perform gate for ",axis.y,  "\n")
  
  if (missing(default_ct_y)) {
    tp_ct_y <-f.roughGate(x, axis.y, thrd = mid.y, tinypeak.removal = initial_tiny_peak, twin.factor = twin.factor.y, alpha = alpha.y)
  } else {
    cat("using default cut for",axis.y,  "\n")
    tp_ct_y <- default_ct_y
  }  
  

    temp_x_pos <- flowDensity(x, c(axis.x, axis.y), position = c(T, NA), gates = c(tp_ct_x, NA))
    temp_x_neg <- flowDensity(x, c(axis.x, axis.y), position = c(F, NA), gates = c(tp_ct_x, NA))
    
    cat("Find cut for ",axis.y, "on", axis.x,  "positive", "\n")
    
    tp_ct_y_right <- f.roughGate(temp_x_pos, axis.y, thrd = mid.y, tinypeak.removal = initial_tiny_peak, twin.factor = twin.factor.y, alpha = alpha.y)
    #tp_ct_y_right <- deGate(temp_x_pos, axis.y, use.upper = T, upper = T, tinypeak.removal = initial_tiny_peak, alpha = alpha.y) # ??
    
    cat("Find cut for ",axis.y, "on", axis.x,  "negative", "\n")
    tp_ct_y_left <- f.roughGate(temp_x_neg, axis.y, thrd = mid.y, tinypeak.removal = initial_tiny_peak, twin.factor = twin.factor.y, alpha = alpha.y)

    if (show.plots == T & temp_x_pos@cell.count > 10) {
      flowPlot(temp_x_pos,
               channels = c(axis.x, axis.y), # xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5),
               title = identifier(x)
      )
    }
    
    if (show.plots == T & temp_x_neg@cell.count > 10) {
      flowPlot(temp_x_neg,
               channels = c(axis.x, axis.y), # xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5),
               title = identifier(x)
      )
    }
    
    
    

    # l.pop <- lapply(children_gates, function(p) {
    #   flowDensity(x, c(axis.x, axis.y), position = p, gates = c(tp_ct_x, tp_ct_y))
    # })
    
    l.pop <- list(
      flowDensity(x, c(axis.x, axis.y), position = children_gates[[1]], gates = c(tp_ct_x, tp_ct_y_right)),
      flowDensity(x, c(axis.x, axis.y), position = children_gates[[2]], gates = c(tp_ct_x, tp_ct_y_left))
    )
    
    
    
    names(l.pop) <- names(children_gates)
    l.gate <- f.generate_gate(l.pop, gate.type = gate.type)
    l.gate_cood <- f.extractGate_boundary(l.gate)
  
  

    
  Gate_colors <- c("blue", "green", "magenta", "brown")
  if (show.plots == T  & nrow(dat) > 10) {
    flowPlot(x,
             channels = c(axis.x, axis.y),  xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5),
             title = identifier(x)
    )
    abline(v = tp_ct_x, h = tp_ct_y, col = "black", lty = 1, lwd = 2)
    if (gate.type == "r" | gate.type == "p") {
      lapply(seq(length(l.gate_cood)), function(x) {
        lines(l.gate_cood[[x]],
              col = Gate_colors[[x]],
              lty = 1, lwd = 2
        )
      })
    }
  }
  
  return(l.gate)
}



# x <- cytoframe_to_flowFrame(fs[[1]]); show.plots = T; parameter = "CD103"; thrd = 2
# only gate on Y axis
f.QC_gating <- function(x, 
                        pop, popname,
                        parameters = c("CD45","CD103"), thrd = 2, alpha = 1/2, position = c(NA,T),
                        show.plots = T){
  
  x <-  cytoframe_to_flowFrame( gs_pop_get_data(x, pop)[[1]])
  
  
  channels <- unlist(lapply(parameters, function(t) {
    channel <- colnames(x)[which(colnames(x) == t)]
    if(length(channel) == 0) {
      channel <- colnames(x)[Find.markers(x, t)]
    }
    
  }))
  
  if(is.na(position[1])) {
    x.rough_cut <-NA
    y.rough_cut <- f.roughGate(x, channels[2], thrd = thrd, alpha = alpha)
  }else{
    x.rough_cut <- f.roughGate(x, channels[2], thrd = thrd, alpha = alpha)   
    y.rough_cut <- f.roughGate(x, channels[2], thrd = thrd, alpha = alpha)    
  }
  
  
  l.pop <- list(
    flowDensity(x, channels, position = position, gates = c(x.rough_cut, y.rough_cut))
  )
  
  if(missing(popname)) {
    popname <- paste0(gsub("_", "", parameters[2]), "+")
  }else{
  }
  names(l.pop) <- popname
  l.gate <- f.generate_gate(l.pop, gate.type = "r")
  l.gate_cood <- f.extractGate_boundary(l.gate)
  
  if (show.plots == T) {
    par(mfrow = c(1, 2))
    flowPlot(x,
             channels = channels, xlim = c(-0.2, 4.5), ylim = c(-0.2, 4.5),
             title = paste0(parameters[2], "by ", parameters[1], "of ", pop),
             title.font.size = 1
    )
    abline(h = y.rough_cut, col = "red")
    abline(h = thrd, col = "grey", lty =2)
    
    flowPlot(x,
             channels = channels, xlim = c(-0.2, 4.5), ylim = c(-0.2, 4.5),
             title = paste(parameters[2], identifier(x)),
             title.font.size = 1
    )
    lapply(seq(length(l.gate_cood)), function(x) {
      lines(l.gate_cood[[x]],
            col =  "blue",
            lty = 1, lwd = 2
      )
    })
    
  }
  
  
  return(l.gate)

}



# x <- fs[[1]]; show.plots = T
f.B_gating <- function(x, show.plots = T) {
  print(identifier(x))

  x.rough_cut <- f.roughGate(x, "B750-A", thrd = 2, alpha = 1 / 100)
  y.rough_cut <- f.roughGate(x, "UV736-A", thrd = 2.5) # print(y.rough_cut)

  if (show.plots == T) {
    par(mfrow = c(1, 4))
    flowPlot(x,
      channels = c("B750-A", "UV736-A"), xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5),
      title = identifier(x)
    )
    abline(h = y.rough_cut, v = x.rough_cut)
  }
  
  CD19pos <- flowDensity(x, c("B750-A", "UV736-A"), position = c(T, NA), gates = c(x.rough_cut, NA))
  
  
  
  pks <- getPeaks(CD19pos, "B750-A")
  if(is.na(pks$Peaks)) {
    cut <- x.rough_cut
  } else {
    if(min(pks$Peaks) < 2){
      if(length(pks$Peaks) == 1) { cut <- deGate(CD19pos, "B750-A", use.upper = T, upper = T)} else {
        cut <- deGate(CD19pos, "B750-A")
        
      }
    } else {
      cut <- x.rough_cut
    }
  }

  
  if(is.infinite(cut)) { cut <- x.rough_cut}

  CD19neg <- flowDensity(x, c("B750-A", "UV736-A"), position = c(F, NA), gates = c(cut, NA))



  y.second_cut <- f.roughGate(CD19neg, "UV736-A", thrd = 2.5, tinypeak.removal = 1 / 100, alpha =1/100)
  # print(y.second_cut)
  cp_NonTB <- flowDensity(x, c("B750-A", "UV736-A"), position = c(F, F), gates = c(cut, y.second_cut))

  # CD19pos <- flowDensity(x, c("B750-A", "UV736-A"), position = c(T,NA), gates = c(x.rough_cut, NA))
  if (show.plots == T) {
    flowPlot(CD19neg, channels = c("B750-A", "UV736-A"), xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5))
    abline(h = y.second_cut)
  }
  

  l.pop <- list(
    cp_NonTB,
    cp_NonTB
  )

  names(l.pop) <- c("NonTB", "B")

  # l.pop_gate <- f.generate_gate(l.pop, gate.type = "r")


  l.pop_gate <- lapply(1:length(l.pop), function(i) {
    g <- apply(l.pop[[i]]@filter, 2, range)
    g[1, ] <- c(-Inf, -Inf) #
    rectangleGate(
      filterId = names(l.pop)[i],
      .gate = g
    )
  })
  names(l.pop_gate) <- names(l.pop)


  l.pop_gate_cood <- f.extractGate_boundary(l.pop_gate)


  if (show.plots == T) {
    flowPlot(x,
      channels = c("B750-A", "UV736-A"), xlim = c(-0.5, 4.5), ylim = c(-0.5, 4.5),
      title = identifier(x)
    )
    lapply(l.pop_gate_cood[1], function(x) {
      lines(x, col = "blue", lty = 1, lwd = 2)
    })
    # abline(v = tp_ct_19, col = "blue")
  }

  return(l.pop_gate)
}





# return a cell identity matrix with row as cell and columns as each category the value is either T or F
f.get_event_identity <- function(gatingset,
                                 color_scheme
                                 # B1_mono_gate = F,
                                 # T_mono_gate = F
) {
  stopifnot(class(color_scheme) == "character")
  obj <- gatingset # gs is required here to extract the population information
  nodes_path <- setNames(gh_get_pop_paths(obj, path = "full"), gh_get_pop_paths(obj, path = "auto")) # extract all the pop names


  cell_identity <- as.data.frame(do.call(cbind, lapply(names(nodes_path), function(x) {
    gh_pop_get_indices(obj, x)
  }))) # extract the cell identity
  colnames(cell_identity) <- gh_get_pop_paths(obj, path = "auto")


  # assign identity based on color_scheme
  coloring_names <- names(color_scheme)
  unique_identity <- cell_identity[, colnames(cell_identity) %in% coloring_names] #subset pops existed in color scheme
  coloring_names <- coloring_names[coloring_names %in% names(unique_identity)]
  unique_identity <- unique_identity[, coloring_names]
  pop <- apply(unique_identity, 1, function(x) {
    colnames(unique_identity)[min(which(x == T))]
  }) # the identity was assigned based on the later column

  pop
}





# improved version
f.prep_plotting_id.list2 <- function(gatingset, id, parent_gate, color_scheme) {
  parent_pop <- gh_pop_get_indices(gatingset, parent_gate)

  id[!parent_pop] <- ""

  id.list <- lapply(names(color_scheme), function(x) {
    which(id == x)
  })

  names(id.list) <- names(color_scheme)

  id.list <- id.list[unlist(lapply(id.list, length)) > 0]
  id.list
}




#' This is a function generate a list of gates based on a list of CellPopulation obj
#' return from flowDensity
#'


# perform gate on a list of cell population
f.generate_gate <- function(l.CellPop, gate.type = "p") {
  if (is.null(names(l.CellPop))) {
    stop("name the populations")
  } else {
    if (gate.type == "p") {
      l.gate <- lapply(1:length(l.CellPop), function(i) {
        polygonGate(
          filterId = names(l.CellPop)[i],
          .gate = l.CellPop[[i]]@filter
        )
      })
    } else if (gate.type == "r") {
      l.gate <- lapply(1:length(l.CellPop), function(i) {
        rectangleGate(
          filterId = names(l.CellPop)[i],
          .gate = apply(l.CellPop[[i]]@filter, 2, range)
        )
      })
    }

    names(l.gate) <- names(l.CellPop)

    return(l.gate)
  }
}

# extract gate from gate list,
#' @@param l.gate a list of gate object
#' for qudrant gate, must get the dp population
#' when dealing with -Inf, assign lower bound to -5 just for visualization
f.extractGate_boundary <- function(l.gate) {
  lapply(l.gate, function(x) {
    gate_type <- paste0(capture.output(x),
      collapse = ""
    )

    if (grepl("Polygonal", gate_type)
    ) {
      x@boundaries
    } else if (grepl("Rectangular", gate_type)
    ) {
      if (grepl("-Inf", gate_type)) {
        if (grepl("NonTB", gate_type)) { # Not the best way
          rbind(
            c(-5, x@max[2]),
            c(x@max[1], x@max[2]),
            c(x@max[1], -5),
            c(4.5, -5),
            c(4.5, 4.5),
            c(-5, 4.5),
            c(-5, x@max[2])
          )
        } else {
          ll <- c(-5, -5)
          lr <- c(x@max[1], -5)
          ur <- c(x@max[1], x@max[2])
          ul <- c(-5, x@max[2])
          rbind(ll, lr, ur, ul, ll)
        }
      } else {
        if (grepl("Inf", gate_type)) {
          x@min
        } else {
          ll <- x@min
          lr <- c(x@max[1], x@min[2])
          ur <- c(x@max[1], x@max[2])
          ul <- c(x@min[1], x@max[2])
          rbind(ll, lr, ur, ul, ll)
        }
      }
    }
  })
}




#' function that handle the gates format
#' it takes in the info passed from csv and part them into a list
# panel_layout = middle_panel_layout
# plotting_gatingset = gs[[1]]
f.plotting_parameters <- function(panel_layout, plotting_gatingset, default_min) {
  
  apply(panel_layout, 1, function(x) {
    gates <- x["gates"]

    if (is.na(x["min.xy"])) {
      if( length(default_min) == 1 ) { min.xy <- c(0, 0) } else {
        min.xy <- c(default_min[ grep(x["x"], default_min$marker_name),1],
                    default_min[ grep(x["y"], default_min$marker_name),1]) 
      }
      
    } else {
      m <- as.character(x["min.xy"])
      min.xy <- as.numeric(strsplit(m, ";")[[1]])
    }

    if (is.na(x["max.xy"])) {
      max.xy <- c(4.5, 4.5)
    } else {
      m <- as.character(x["max.xy"])
      max.xy <- as.numeric(strsplit(m, ";")[[1]])
    }




    if (!is.na(gates)) {
      gates <- gsub(" ", "", gates)
      gates <- strsplit(gates, ";")[[1]]
      
      l.gates <- lapply(gates, function(g) {
        exists_gates <-gh_get_pop_paths(plotting_gatingset, path="auto")
        if(g %in% exists_gates) {
          gate <- gh_pop_get_gate(plotting_gatingset, g)
        } else {
          cat("No gate was found\n")
          ""
        }
      }) # get gate object

      l.gates <- f.extractGate_boundary(l.gates)
    } else {
      (l.gates <- "")
    }

    
    if(!is.na(x["highlight"][1])) {
      hl <- gsub(" ", "", x["highlight"])
      hl <- strsplit(hl, ";")[[1]]
    }else{
      hl <-NA
    }
    

    list(
      axis = c(x["x"], x["y"]),
      Parent_gate = x["parent"],
      type = x["type"],
      text = x["text"],
      gates = l.gates,
      min.xy = min.xy,
      max.xy = max.xy,
      z = x["z"],
      highlight = hl,
      subtitle = x["subtitle"],
      subtitle.color = x["subtitle.color"]
    )
  })
}




f.nearest_to_reference <- function(cuts, reference) {
  cuts[which.min(abs(cuts - reference))]
}



f.add_w_to_cf <- function(cf, transform_path ) {
  para <- read.csv(transform_path, row.names = 1)
  p <-parameters(cf)
  p@data$trans <-NA
  
  p@data$trans <- para$w[match(p@data$name, rownames(para))]
  p@varMetadata <- rbind(p@varMetadata, trans = "Transformation Parameter")
  parameters(cf) <- p
  
  ff <-cytoframe_to_flowFrame(cf)
  ff@parameters <- p
  
  ff
  
  }


library(scales)
#' Gate Beads in Flow Cytometry Data
#'
#' This function performs gating on bead populations in flow cytometry data. It identifies
#' bead populations based on forward scatter (FSC-A) and side scatter (SSC-A) channels,
#' applies gates to isolate bead populations, and optionally visualizes the results.
#'
#' @param x A flow cytometry data object (e.g., a `flowFrame` or similar object).
#' @param FCS.max Numeric. The maximum value for the FSC-A channel. Values above this
#'                threshold are truncated. Default is `500000`.
#' @param SSC.max Numeric. The maximum value for the SSC-A channel. Values above this
#'                threshold are truncated. Default is `400000`.
#' @param show.plots Logical. If `TRUE`, plots are generated to visualize the gating
#'                   process and results. Default is `TRUE`.
#'
#' @return A list containing:
#' \itemize{
#'   \item `beads`: A `flowFrame` object representing the gated bead population.
#'   \item `beads_plus`: A `flowFrame` object representing the gated bead population
#'                       with additional filtering.
#'   \item `l.gates`: A list of gates generated for the bead populations.
#' }
#'
#' @details
#' The function performs the following steps:
#' 1. Truncates FSC-A and SSC-A values to the specified maximum thresholds.
#' 2. Identifies bead populations using the `flowDensity` package.
#' 3. Applies additional gating to refine the bead populations.
#' 4. Optionally generates plots to visualize the gating process.
#' 5. Returns the gated bead populations and associated gates.
#'
#' @examples
#' \dontrun{
#' # Load a flow cytometry dataset
#' library(flowCore)
#' data <- read.FCS("path_to_fcs_file.fcs")
#'
#' # Gate bead populations
#' bead_gates <- f.Gate_beads(data, FCS.max = 500000, SSC.max = 400000, show.plots = TRUE)
#'
#' # Access gated bead populations
#' beads <- bead_gates$beads
#' beads_plus <- bead_gates$beads_plus
#' gates <- bead_gates$l.gates
#' }
#'
#' @export
#' 
f.Gate_beads <- function(x, FCS.max =500000, SSC.max = 400000, show.plots = T ) {
  
  #exprs(x)[, "FSC-A"] <- pmin(FCS.max, exprs(x)[, "FSC-A"])
  #exprs(x)[, "SSC-A"] <- pmin(SSC.max, exprs(x)[, "SSC-A"])
  print(identifier(x))
  t <- gsub(".fcs", "", identifier(x))
  t <- gsub("Compensation Controls_", "", t)
  x <- flowDensity(x, c("FSC-A", "SSC-A"), 
                   position = c(F, F), 
                   gates = c(FCS.max, SSC.max))
  x <- flowDensity(x, c("FSC-A", "SSC-A"), 
                   position = c(T, T), 
                   gates = c(10000, 50000))
  
  
  
  
  if(show.plots == T) {
    par(mfrow = c(1,5))
    flowPlot(x, channels = c("FSC-A", "SSC-A"), title = t)
    
  }
  tinypeak <- 1/2
  ssc_peaks <- getPeaks(x, "SSC-A", tinypeak.removal = tinypeak) # make sure to get major peaks
  fsc_peaks <- getPeaks(x, "FSC-A", tinypeak.removal = tinypeak) # make sure to get major peaks
  
  
  if(length(ssc_peaks$Peaks) >1) {
    ssc_cut_1 <- deGate(x, "SSC-A", tinypeak.removal = tinypeak)
  }else {
    ssc_cut_1 <- deGate(x, "SSC-A", tinypeak.removal = tinypeak, use.upper =  T, upper = T)
  }
  
  if(length(fsc_peaks$Peaks) >1) {
    fsc_cut_1 <- deGate(x, "FSC-A", tinypeak.removal = tinypeak)
  }else {
    fsc_cut_1 <- deGate(x, "FSC-A", tinypeak.removal = tinypeak, use.upper =  T, upper = T)
  }
  
  beads <- flowDensity(x, c("FSC-A", "SSC-A"), 
                       position = c(F, F), 
                       gates = c(fsc_cut_1, ssc_cut_1))
  
  beads_plus <- flowDensity(x, c("FSC-A", "SSC-A"), 
                            position = c(T, T), 
                            gates = c(fsc_cut_1, ssc_cut_1))
  
  
  beads_cuts_2 <- lapply(c(F,T), function(l) {
    c( `SSC` = deGate(beads, "SSC-A", tinypeak.removal = tinypeak, use.upper =  T, upper = l),
       `FSC` =deGate(beads, "FSC-A", tinypeak.removal = tinypeak, use.upper =  T, upper = l))
    
  })
  
  beads_plus_cuts_2 <- lapply(c(F,T), function(l) {
    c( `SSC` = deGate(beads_plus, "SSC-A", tinypeak.removal = tinypeak, use.upper =  T, upper = l),
       `FSC` =deGate(beads_plus, "FSC-A", tinypeak.removal = tinypeak, use.upper =  T, upper = l))
    
  })
  
  
  if(show.plots == T) {
    flowPlot(beads, channels = c("FSC-A", "SSC-A"))
    abline(h = beads_cuts_2[[1]][1], v =  beads_cuts_2[[1]][2], col = "red")
    abline(h = beads_cuts_2[[2]][1], v =  beads_cuts_2[[2]][2], col = "red")
    flowPlot(beads_plus, channels = c("FSC-A", "SSC-A"))
    abline(h = beads_plus_cuts_2[[1]][1], v =  beads_plus_cuts_2[[1]][2], col = "red")
    abline(h = beads_plus_cuts_2[[2]][1], v =  beads_plus_cuts_2[[2]][2], col = "red")
  }
  
  beads <-    flowDensity(beads, c("SSC-A", "FSC-A"), position = c(T, T), gates = beads_cuts_2[[1]] )
  beads <-    flowDensity(beads, c("SSC-A", "FSC-A"), position = c(F, F), gates = beads_cuts_2[[2]] )
  
  beads_plus <-    flowDensity(beads_plus, c("SSC-A", "FSC-A"), position = c(T, T), gates = beads_plus_cuts_2[[1]] )
  beads_plus <-    flowDensity(beads_plus, c("SSC-A", "FSC-A"), position = c(F, F), gates = beads_plus_cuts_2[[2]] )
  
  
  if(show.plots == T) {
    flowPlot(beads, channels = c("FSC-A", "SSC-A"))
    flowPlot(beads_plus, channels = c("FSC-A", "SSC-A"))
    
  }
  l.pops <- list(`beads` = beads,  `beads_plus`= beads_plus)
  l.gates <- f.generate_gate(l.pops, gate.type = "r")
  
  
  return(l.gates)
  
}







f.inflection_cuts <- function(x, channel) {
  
  for (i in c(10,8,6,4,2,1)) {
    p<- getPeaks(x, channel, tinypeak.removal = 1/i)
    if(length(p$Peaks) ==1) {break}
  }
  
  c(
    deGate(x, channel, tinypeak.removal = 1/i, use.upper = T, upper = F),
    deGate(x, channel, tinypeak.removal = 1/i, use.upper = T, upper = T)
  )
} 

# 
# i=1
# x <- cytoframe_to_flowFrame(fs[[i]])
# channel <- channels[i]
# show.plots <- T
# channel.max <- 262143
# xlim = c(-2,8)

f.Gate_single_stained_control <- function(x,
                                          channel, channel.max = 262143 , 
                                          threshold = 0.5,
                                          show.plots = T,
                                          xlim = c(-2,8),
                                          pop_size_limit = 5000
                                          ) {
  t <- gsub(".fcs", "", identifier(x))
  t <- gsub("Compensation Controls_", "", t)
  
  
  x_frame <- cytoframe_to_flowFrame(x)
  
  
  #channel <- colnames(x)[grep(channel,colnames(x))]
  exprs(x_frame)[, channel] <- asinh(exprs(x_frame)[, channel] / 400)
  channel.max_asinh <- asinh(channel.max / 400)
  
  
  # x <- flowDensity(x_frame, c(channel, "SSC-A"),
  #                  position = c(F, NA),
  #                  gates = c(channel.max_asinh, NA)
  # )
  
  x <- x_frame
  cuts <- f.inflection_cuts(x, channel)
  x <- flowDensity(x, c(channel, "SSC-A"), position = c(T, NA), gates = c(cuts[1], NA))
  x <- flowDensity(x, c(channel, "SSC-A"), position = c(F, NA), gates = c(cuts[2], NA))
  
  peaks <- getPeaks(x, channel, tinypeak.removal = 1 / 10, twin.factor = 1)
  
  if (length(peaks$Peaks) == 1) { # only one peaks   
    
    cuts <- f.inflection_cuts(x, channel)
    cuts <-pmin(cuts, channel.max_asinh)
    
    if (peaks$Peaks > threshold) { #only one positive peak
      
      pos <- flowDensity(x, c(channel, "SSC-A"), position = c(T, NA), gates = c(cuts[1], NA))
      pos <- flowDensity(pos, c(channel, "SSC-A"), position = c(F, NA), gates = c(cuts[2], NA))

      neg <- flowDensity(x, c(channel, "SSC-A"), position = c(F, NA), gates = c(-Inf, NA))
      
      if (show.plots == T) {
        par(mfrow = c(1, 4))
        flowPlot(x_frame, channels = c(channel, "SSC-A"), title = t, xlim = xlim)

        #abline(v = cuts, col = "black")
        #abline(v = cuts_2_neg, col = "blue")
        abline(v = cuts, col = "red")
        abline(v = channel.max_asinh, col = "blue", lw=2)
        flowPlot(neg, channels = c(channel, "SSC-A"), title = t, xlim = xlim)
        abline(v = channel.max_asinh, col = "blue", lw=2)
        
        flowPlot(pos, channels = c(channel, "SSC-A"), title = t, xlim = xlim)
        abline(v = channel.max_asinh, col = "blue", lw=2)
        
      }
      
      
    } else {
      pos <- flowDensity(x, c(channel, "SSC-A"), position = c(F, NA), gates = c(-Inf, NA))
      #pos$filter <- cbind(c(-Inf, -Inf,-Inf,-Inf), c(-Inf, -Inf,-Inf,-Inf))
      
      neg <- flowDensity(x, c(channel, "SSC-A"), position = c(T, NA), gates = c(cuts[1], NA))
      neg <- flowDensity(neg, c(channel, "SSC-A"), position = c(F, NA), gates = c(cuts[2], NA))
      
      
      if (show.plots == T) {
        par(mfrow = c(1, 4))
        flowPlot(x_frame, channels = c(channel, "SSC-A"), title = t, xlim = xlim)

        abline(v = cuts, col = "blue")
        abline(v = channel.max_asinh, col = "green")
      }
      
      
    }
  } else { length(peaks$Peaks) >=2 
      
    
    if(max(peaks$Peaks) <= threshold) { #no positive peaks
      cuts <- cuts_2_pos <- NULL
      cuts_2_neg <- f.inflection_cuts(x, channel)
      
      
      cuts_between <- deGate(x, channel, twin.factor =1)
      
      
      #neg <- flowDensity(x, c(channel, "SSC-A"), position = c(T, NA), gates = c(cuts_2_neg[1], NA))
      neg <- flowDensity(x, c(channel, "SSC-A"), position = c(T, NA), gates = c(cuts_between, NA))
      neg <- flowDensity(neg, c(channel, "SSC-A"), position = c(F, NA), gates = c(cuts_2_neg[2], NA))
      pos <- flowDensity(x, c(channel, "SSC-A"), position = c(F, NA), gates = c(-Inf, NA))
      
      
      
    } else if(min(peaks$Peaks) >= threshold ) { # all positive peaks
      
      cuts <- cuts_2_neg <- NULL
      cuts_2_pos <- f.inflection_cuts(x, channel)
      cuts_2_pos <-pmin(cuts_2_pos, channel.max_asinh)
      cuts_between <- deGate(x, channel, twin.factor =1)
      
      #pos <- flowDensity(x, c(channel, "SSC-A"), position = c(T, NA), gates = c(cuts_2_pos[1], NA))
      pos <- flowDensity(x, c(channel, "SSC-A"), position = c(T, NA), gates = c(cuts_between, NA))
      pos <- flowDensity(pos, c(channel, "SSC-A"), position = c(F, NA), gates = c(cuts_2_pos[2], NA))
      
      if(pos@cell.count < pop_size_limit) {
        pos <- flowDensity(x, c(channel, "SSC-A"), position = c(T, NA), gates = c(cuts_2_pos[1], NA))
        pos <- flowDensity(pos, c(channel, "SSC-A"), position = c(F, NA), gates = c(cuts_2_pos[2], NA))
      }
      
      neg <- flowDensity(x, c(channel, "SSC-A"), position = c(F, NA), gates = c(-Inf, NA))
      
      
        }  else {
      cuts <- deGate(x, channel, tinypeak.removal = 1 / 10, twin.factor = 1)
      
      
      # rough cut remove
      neg <- flowDensity(x, c(channel, "SSC-A"), position = c(F, NA), gates = c(cuts, NA))
      pos <- flowDensity(x, c(channel, "SSC-A"), position = c(T, NA), gates = c(cuts, NA))
      
      cuts_2_neg <- f.inflection_cuts(neg, channel)
      
      # fine cuts
      neg <- flowDensity(neg, c(channel, "SSC-A"), position = c(T, NA), gates = c(cuts_2_neg[1], NA))
      neg <- flowDensity(neg, c(channel, "SSC-A"), position = c(F, NA), gates = c(cuts_2_neg[2], NA))
      
      
      cuts_2_pos <- f.inflection_cuts(pos, channel)
      cuts_2_pos <-pmin(cuts_2_pos, channel.max_asinh)
      
      
      pos <- flowDensity(pos, c(channel, "SSC-A"), position = c(T, NA), gates = c(cuts_2_pos[1], NA))
      pos <- flowDensity(pos, c(channel, "SSC-A"), position = c(F, NA), gates = c(cuts_2_pos[2], NA))
      
    }
    
    
    
    
    if (show.plots == T) {
      par(mfrow = c(1, 4))
      flowPlot(x_frame, channels = c(channel, "SSC-A"), title = t, xlim = xlim)
      abline(v = cuts, col = "black")
      abline(v = cuts_2_neg, col = "blue")
      abline(v = cuts_2_pos, col = "red")
      abline(v = cuts_between, col = "yellow")
      abline(v = channel.max_asinh, col = "blue", lw=2)
      flowPlot(neg, channels = c(channel, "SSC-A"), title = t, xlim = xlim)
      abline(v = channel.max_asinh, col = "blue", lw=2)
      flowPlot(pos, channels = c(channel, "SSC-A"), title = t, xlim = xlim)
      abline(v = channel.max_asinh, col = "blue", lw=2)
      
      
    }
    
    
  }
  
  
  l.pops <- list(`pos` = pos,  `neg`= neg)
  l.gates <- f.generate_gate(l.pops, gate.type = "r")
  
  l.gates <- lapply(l.gates, function(x) {
       mn <- x@min[channel]
       mx <- x@max[channel]
    
      if(!is.infinite(mn)) {x@min[channel] <- sinh(mn) *400}
       if(!is.infinite(mx)) {x@max[channel] <- sinh(mx) *400}
     x
  })
  
  
  return(l.gates)
  

}






#pops <- c("T; NK; B")
#colorkey <- autoflow_colors_CLT$autoflow_30c_T
f.plot_color_legend <- function(pops, colorkey) {
  pops <- gsub(" ", "", unlist(strsplit(pops, ";")))
  cls <- rev(colorkey[names(colorkey) %in% pops])
  
  
  
  n <-9
  legs <- rep("white", (n-length(cls)))
  names(legs) <- rep("", (n-length(cls)))
  
  legs <- c(cls, legs)
  plot(1, 1, xlim = c(1, 5.5), ylim = c(0, n), type = "n", ann = FALSE, axes = FALSE)
  for(i in 1:(n-1)) {
    points(1.1, n - i, cex = 3, col = legs[i], pch = 16)
    text(1.5, n - i, cex = 1.5, names(legs) [i], adj = 0)
  }
  
  
}
  
  
  









#' This is a helper function for plotting a plots panel include multiple plots
#' @param obj gatingset
#' @layout layout  description, its a dataframe ol

# obj = gs[[1]]
# layout = middle_panel_layout
# transform_csv_path = transformation.parameters.path
# colors = autoflow_colors$autoflow_30c_T
# width = 12
# height = 12
# postfix = "middle"
#colors = autoflow_colors$autoflow_30c_B
#res = NA

f.plot_panel <- function(obj, DR_flowframe, layout, width = 12, height = 12, # nrow = 5, ncol = 5,
                         rendering.type = "nonCairo", res = NA, Gatefinder = NA,
                         show.time = F, save.png = F, path.png = "", 
                         transform_csv_path = transformation.parameters.path, 
                         colors, 
                         gates.color = "grey50",
                         global.min, # = autoflow_colors$autoflow_30c_T, highlight,
                         postfix = NULL) {
  
  
  plotting_gs <- obj
  if (missing(DR_flowframe)) {
    dat <- gs_pop_get_data(plotting_gs, "root")[[1]]
  } else {
    dat <- DR_flowframe
  }

  # add transformation info to the flowframe
  dat <- f.add_w_to_cf(dat, transform_path = transform_csv_path) 

  print(identifier(dat))
  id_vec <- f.get_event_identity(plotting_gs, color_scheme = colors)
  if(missing(global.min)){global.min <- NA}
  plotting_parameters <- f.plotting_parameters(layout, plotting_gs, global.min)
  sample_name <- gsub(".fcs", "", sampleNames(plotting_gs))
  



  nrow <- max(layout$row)
  ncol <- max(layout$colum)

  if("annotation" %in% colnames(layout)) {
    annotations <- layout$annotation
  } else {
    annotations <- rep(NA, nrow(layout))
  }
  
  
  
  
  
  fig.width <- width
  fig.height <- height

  dev.new(#width = fig.width, height = fig.height, 
          noRStudioGD = T)
  # path.png <- "/Users/xuehai/Downloads/"

  if (save.png == T) {
    output_filenames <- paste0(path.png, "/", sample_name, postfix, ".png")
    if (rendering.type == "Cairo") {
      Cairo::CairoPNG(
        filename = output_filenames,
        width = fig.width, height = fig.height, 
        units = "in", res = 300
      )
    } else {
      png(
        filename = output_filenames,
        width = fig.width, height = fig.height, 
        units = "in", res = 300
      )
    }
  }



  par(
    mar = c(3.5, 3.5, 1.5, 2), # bottom, left side.....
    mgp = c(2, 0.5, 0), # axis label, tick labels, better not move
    oma = c(0, 2, 0, 0),
    mfrow = c(nrow, ncol) # ,
    # xpd = NA
  )
  dev.control("enable") # this is the key to make sure the plots were saved

  start_time <- Sys.time()



 
 
 if ( (!missing(DR_flowframe)) & sum(grepl("UMAP|PC", colnames(dat)))==0 ) { #  
 
 }else{ 
   cat("[Plotting]", "\n")
   for (i in 1:length(plotting_parameters)
   ) {
     if (plotting_parameters[[i]]$type == "Text") {
       plot(1, 1, xlim = c(0, 2), ylim = c(0, 2), type = "n", ann = FALSE, axes = FALSE)
       text(0, 1,
            cex = 2,
            plotting_parameters[[i]]$text,
            adj = 0
       )
     } else if (plotting_parameters[[i]]$type == "legend") {
       #colorkey <- autoflow_colors_CLT$autoflow_30c_T
       f.plot_color_legend(plotting_parameters[[i]]$text, colorkey = colors)
     } else if (plotting_parameters[[i]]$type == "GateFinder") {
       
       # 
       # plot(1, 1, xlim = c(0, 2), ylim = c(0, 2), type = "n", ann = FALSE, axes = FALSE)
       # text(0, 1,
       #      cex = 2,
       #      plotting_parameters[[i]]$text,
       #      adj = 0
       # )
       plot(Gatefinder[[plotting_parameters[[i]]$text]])
       
       
     } else {
       # highlight B subsets
       plotting_id_list <- f.prep_plotting_id.list2(
         gatingset = plotting_gs,
         id = id_vec,
         parent_gate = plotting_parameters[[i]]$Parent_gate,
         color_scheme = colors
       )
       
       tt <- strsplit(plotting_parameters[[i]]$Parent_gate, "/")[[1]]
       tt <- tt[length(tt)]
       
       #tt <- StrExtract(plotting_parameters[[i]]$Parent_gate, "/", 2)
       print(plotting_parameters[[i]]$min.xy)
       PlotFlowFrame(dat,
                     plotting_parameters[[i]]$axis,
                     pop_idex.list = plotting_id_list,
                     colormap_marker = plotting_parameters[[i]]$z,
                     population_colors = colors,
                     title = tt, #plotting_parameters[[i]]$Parent_gate,
                     title.font.size = 1.8,
                     subtitle = plotting_parameters[[i]]$subtitle,
                     subtitle.color = plotting_parameters[[i]]$subtitle.color,
                     min.xy = plotting_parameters[[i]]$min.xy,
                     max.xy = plotting_parameters[[i]]$max.xy,
                     log_ssc = T,
                     res = res,
                     #dens_res = 128,
                     highlight_pops = plotting_parameters[[i]]$highlight,
                     gates = plotting_parameters[[i]]$gates,
                     gates_color = gates.color,
                     plotting_background = F, interpolate = F,
                     rare_on_top = T, auto_assign_color = F,
                     rgba_blending = T, old_method = F,  #shownumbers = T,
                     fulllength_label = F, label_size =1.2
       )
       
       #print(i)
       if(!is.na(annotations[i])) {
         ann <- strsplit(annotations[i], ";")[[1]]
         ann <- gsub(" ", "", ann)
         f.add_pct_to_plot (
           f.percent_annotation(ann, denominator =  layout$parent[i], gh = plotting_gs)  
           
         )
       }
       
       
       
       
       
       
     }
   }
 }
 
 
 


  end_time <- Sys.time()
  if (show.time == T) {
    print(end_time - start_time)
  }

  if (save.png == T) {
    dev.off()
  } else {
    end_time <- Sys.time()

    pl <- recordPlot()
    dev.off()
    return(pl)
  }
}



# inputfolder <- "/Users/flowuser/Documents/Debug/FCS/test"
# a folder path as input
#' out put list with a tables and a vector of the filenames of fcs files need to be analyzed

f.get_meta_data <- function(inputfolder) {
  filenames <- list.files(inputfolder, pattern = ".fcs", full.names = F) # dir()

  if (length(filenames) == 0) {
    filenames <- Data_Table.MOVE <- NULL
  } else {
    # positions <- StrExtract(StrExtract(filenames, ".fcs", 1), "_", 4)
    # positions[which(nchar(positions) != 3)] <- "TM" # Assign tube mode in case there is no well poistion
    tube <- StrExtract(filenames, "_", 3)
    tube <- gsub(".fcs", "", tube) # if no well poistion is there, ".fcs" will be included which need to be removed

    Data_Table <- data.frame(
      BF = StrExtract(filenames, " ", 1),
      PATIENT = paste(StrExtract(StrExtract(filenames, " ", 3), "_", 1), StrExtract(StrExtract(filenames, " ", 3), "_", 2), sep = "_"),
      TISSUE = StrExtract(filenames, " ", 2),
      TUBE = tube, # StrExtract(filenames, "_",3),
      # POSITION = positions, # StrExtract(StrExtract(filenames, ".fcs",1), "_", 4),
      FCS = paste0(inputfolder, "/", filenames),
      FCS_NAME = filenames,
      # SIZE = file.size(filenames),
      ANALYZE = F, # Binary, will be assined TRUE if all 3 tubes were found for analysis
      NOTE = " ",
      stringsAsFactors = F
    )

    l.Data_Table <- split(Data_Table, Data_Table$BF)


    # check the fcs file per BF number
    l.Data_Table <- lapply(l.Data_Table, function(x) {
      b.index <- which(grepl("tubeB", x$TUBE, ignore.case = T))
      t.index <- which(grepl("tubeT", x$TUBE, ignore.case = T))



      if (length(b.index) >= 1 & length(t.index) >= 1) { # if ther is
        x[c(b.index[1], t.index[1]), ]$ANALYZE <- T
        x[c(b.index[1], t.index[1]), ]$NOTE <- "2 tubes"
      } else if (length(b.index) >= 1 & length(t.index) == 0) {
        x[c(b.index[1]), ]$ANALYZE <- T
        x[c(b.index[1]), ]$NOTE <- "Single B"
      } else if (length(b.index) == 0 & length(t.index) >= 1) {
        x[c(t.index[1]), ]$ANALYZE <- T
        x[c(t.index[1]), ]$NOTE <- "Single T"
      }

      return(x)
    })

    Data_Table <- do.call(rbind, l.Data_Table)

    Data_Table$`CELL NUMBER` <- as.numeric(
      unlist(
        lapply(Data_Table$FCS, function(x) {
          # read.FCS(x)@description$`$TOT`
          read.FCSheader(x, keyword = "$TOT")
        })
      )
    )

    Data_Table$`BCCA` <- as.numeric(
      unlist(
        lapply(Data_Table$FCS, function(x) {
          # read.FCS(x)@description$`$TOT`
          read.FCSheader(x, keyword = "BCCA#")
        })
      )
    )


    # Data_Table$'CELL NUMBER' <- format(Data_Table$'CELL#', nsmall=0, big.mark=",")

    Data_Table <- Data_Table[, c(1, 2, 10, 3, 4, 6, 9, 8, 7)]


    return(
      Data_Table
      # list(
      # display_table = Data_Table,
      # files_tobe_analzyed = Data_Table$`FCS NAME`[Data_Table$ANALYZE]
      # )
    )
  }
}








#' this the main function that take all the setting then perform the automated gating
#'
#'
#' This version do not require SETTINGS.R.
#' The original version need to sourcing the SETTINGS.R file multiple times due to each piece of code reset the environemnt.
#' removing the SETTINGS.R because some hicups in the GUI
#'
#' This function will copy and paste the FCS files from the raw foder into analyze foder
#' will also copy paste and update the scripts used for the automated analysis
#'
#' apply modifications to the script to makke it case specfic
#'
#' Perform gating, genrate result html and move the html to dedicated folder
#'
#' @param AUTOFLOW_destination: destination fo the AUTOFLOW program
#' @param raw_fcs_folder original folders contain fcs files needed to be analzyed description
#' @param
#' @param path.output: string object; directory to the folders storing the all the temporary files.
#' @param batch: string object; usually the current date, it is needed to create a run-specific folder "Results_folder"
#' under the path.output to store temporary files for each run
#' @param max_cell_number_to_plot: Maximun cell number to be plotted
#' @param keeptimepfiles: 
#' @param max_analyzing_number maximum cell number to analyze
#' @param Data_Table.MOVE a table contains files moved from raw data file
#' @param filenames filenames of files to be analyzed
#' @param T_function Function used to analysis T tube
#' @param REPORT_function function used to generate html report
#' @param DPI dpi, defaul 72 (was 36 for 14c data)

#' Intermediate parameters
#' @param temp_fcs_folder: folder under "Results_folder" for temporarily storing FCS files

# AUTOFLOW_destination = paste0("~/Developer/AUTOFLOW32", "/")
# path.output = "/Users/flowuser/Documents/Reports_OUT"
# batch = "202502"
# raw_fcs_folder = "/Users/flowuser/Documents/Debug/FCS/test"
# max_analyzing_number = 5000000
# Data_Table <- f.get_meta_data(raw_fcs_folder)
# B_function = "32C_B_AUTOGATING_FD_2025A5.Rmd"
# T_function = "32C_T_AUTOGATING_FD_2025A5.Rmd"
# REPORT_function = "BCCA_LPD_AUTOGATING_32.Rmd"
# report_only_folder = NULL
# DPI = 72

f.autoflow3 <- function(AUTOFLOW_destination,
                        raw_fcs_folder,
                        Data_Table,
                        max_analyzing_number,
                        keeptimepfiles,
                        B_function,# = "32C_B_AUTOGATING_FD_2025A5_CLF.Rmd",
                        T_function,# = "32C_T_AUTOGATING_FD_2025A5.Rmd",
                        B_comp,
                        T_comp,
                        REPORT_function = "BCCA_LPD_AUTOGATING_32.Rmd",
                        report_only_folder = NULL,
                        stats_folder = NULL,
                        path.output,
                        batch,
                        DPI = 72) {
  
  Results_folder <- paste(path.output, batch, sep = "/") # creating run-specific folder
  temp_scripts_folder <- paste(Results_folder, "scripts", sep = "/") # Folder under Results_folder for store temporary scripts
  temp_fcs_folder <- paste(Results_folder, "fcs", sep = "/") # folder with FCS files coppying from the raw fcs files


  # Creating folder for keeping the automated pipeline outputs
  suppressWarnings(dir.create(temp_scripts_folder, recursive = T)) #
  suppressWarnings(dir.create(temp_fcs_folder, recursive = T))
  suppressWarnings(dir.create(paste(temp_scripts_folder, "reports_scripts", sep = "/"), recursive = T))

  suppressWarnings(dir.create(paste(Results_folder, "gating", sep = "/"), recursive = T))
  suppressWarnings(dir.create(paste(Results_folder, "stats/counts", sep = "/"), recursive = T))
  suppressWarnings(dir.create(paste(Results_folder, "workspace", sep = "/"), recursive = T))
  suppressWarnings(dir.create(paste(Results_folder, "data", sep = "/"), recursive = T))
  suppressWarnings(dir.create(paste(Results_folder, "reports", sep = "/"), recursive = T))
  suppressWarnings(dir.create(paste(Results_folder, "html", sep = "/"), recursive = T))
  suppressWarnings(dir.create(paste(Results_folder, "png", sep = "/"), recursive = T))
  suppressWarnings(dir.create(paste(Results_folder, "comp", sep = "/"), recursive = T))
  
  ### display the folder path ========================
  message(paste0("RAW fcs folder: ", raw_fcs_folder, "\n"))
  message(paste0("Result folder: ", Results_folder, "\n"))

  # Copy/Move the FCS files ========================
  Copy <- T

  ## remove the existing files under the temp_fcs_folder =============
  if (length(dir(paste0(temp_fcs_folder))) > 0) {
    file.remove(list.files(paste0(temp_fcs_folder), pattern = ".fcs", full.names = T))
  }

  filenames <- Data_Table[Data_Table$ANALYZE == T, "FCS_NAME"]
  message("\nFollowing fcs files will be copied and analyzed:\n")
  message(paste(filenames, collapse = "\n"))

  if (Copy == T) {
    message("Copy FCS files!!\n")

    file.copy(
      from = paste0(raw_fcs_folder, "/", filenames),
      to = paste0(temp_fcs_folder), overwrite = F
    )
  } else {
    filesstrings::file.move(
      files = paste0(raw_fcs_folder, "/", filenames),
      destinations = paste0(temp_fcs_folder), overwrite = F
    )
  }
  
  #B_comp_path <- B_comp$PATH
  B_comp_filename <- B_comp$FILE
  B_comp_path_local <- paste(Results_folder, "comp", B_comp_filename, sep = "/")
  
  #T_comp_path <- T_comp$PATH
  T_comp_filename <- T_comp$FILE
  T_comp_path_local <- paste(Results_folder, "comp", T_comp_filename, sep = "/")

  ## Copy comp files
  if (Copy == T) {
    message("Copy COMPENSTAION files!!\n")
    
    file.copy(
      from = B_comp$PATH,
      to = B_comp_path_local
    )
    
    file.copy(
      from = T_comp$PATH,
      to = T_comp_path_local
    )
    
    
    
  } else {
    filesstrings::file.move(
      files = paste0(raw_fcs_folder, "/", filenames),
      destinations = paste0(temp_fcs_folder), overwrite = F
    )
  }
  
  
  # read and parse the fcs files ========================
  # select files to analzye according to data table
  # files <- dir(temp_fcs_folder)

  # Create a file list for the algorithm to cycle through the analysis
  Data_Table_T <- Data_Table[Data_Table$ANALYZE == T, ]
  Data_Table_T$NAME <- paste(Data_Table_T$BF, Data_Table_T$TISSUE, Data_Table_T$PATIENT, Data_Table_T$NOTE, sep = " ")
  Data_Table_T$NAME <- gsub(" ", "_", Data_Table_T$NAME)
  file.list <- split(Data_Table_T, Data_Table_T$NAME) # one patient per element in the list


  # Analyse the sample by sample  ===============
  ptm <- proc.time()

  Report_filenames <- NULL

  for (s in names(file.list)) {
  log <- tryCatch(
      {
        case <- file.list[[s]]

        l.case.data <- list() # to collect the result data path for each row, its possible that there are more than two B or T (repeated cases)

        for (cs in 1:nrow(case)) { # changed cs to i so it will not be affected by B and T scriptsi
          
          cat("Starting analyze", case$TUBE[cs], "for case:", s, "...\n")

          if (case$TUBE[cs] == "tubeB") {
            
            tryCatch(
              {
            
            cat("Processing tube B fcs file for sample: ", case$BF[cs], "...\n")
            #print(case$TUBE[cs])


            ### Update B1 script content ===============
            # condition to select normal or QC analyzing file
            if(grepl("QC", s, ignore.case = T) & grepl("CHEX", s, ignore.case = T) ) {
              
              warning("Loading QC Analyzing template")
              file_contents <- readLines(paste0(AUTOFLOW_destination, "/R/", 
                                                gsub(".Rmd", "", B_function),
                                                "_QC.Rmd"))
            } else {
              file_contents <- readLines(paste0(AUTOFLOW_destination, "/R/", 
                                                B_function))
              }
                
                
                
            
            
            file_contents <-  file_contents[1: (which(file_contents == "# DEV PART #####################################################################" )-1)]
            
            fcs.names <- case$FCS_NAME[cs]

            
            ## a list of content that need to be updated
            list.content_update <- list (
              `afpath`= c("~/Developer/AUTOFLOW32/", AUTOFLOW_destination),
              `filename` = c("XXX_FILENAME_B_XXX", fcs.names),
              `filepath` = c("XXX_FILEPATH_B_XXX", paste(temp_fcs_folder, fcs.names, sep = "/")),
              `comppath` = c("XXX_COMP_PATH_B_XXX", B_comp_path_local),
              `compfile` = c("XXX_COMP_FILE_B_XXX", B_comp_filename),
              #`pngfolder` = c("XXX_PNG.PATH_XXX", paste(Results_folder, "png", sep = "/")),
              `path.output` = c("/Users/flowuser/Documents/Debug/OUTPUT", path.output),
              `batchfolder` = c("XXX_BATCHFOLDER_XXX", batch),
              `max_analyzing_number` = c("85678935", max_analyzing_number)
            )
            
            ## Update B script content ===============
            for (l in seq(length(list.content_update))) {
              file_contents <- gsub(
                x = file_contents, 
                pattern = list.content_update[[l]][1],
                replacement = list.content_update[[l]][2]
              )
            }

            ### save the updated T analysis file
            temp_t_script_path <- paste(temp_scripts_folder, paste(batch, case$BF[cs], "B.Rmd", sep = "_"), sep = "/")
            cat(file_contents, file = temp_t_script_path, sep = "\n")
            
            # Run B script ==================
            rmarkdown::render(temp_t_script_path, #output_file = "aa.html",
                              output_dir = paste(Results_folder, "html", sep = "/"), envir = new.env(),
                              clean = T)
            


          }, error = function(e) {
            s <- paste0("FAILED_", s, "TubeB.txt")
            sink(paste0(paste(Results_folder, "reports", sep = "/"), "/", "FAILED_", s, ".txt"))
            cat("ERROR :", conditionMessage(e), "\n")
            sink()
          }
         )
            data_path <- paste(Results_folder, "data", paste0("l.B.data.", StrExtract(fcs.names, ".fcs", 1), ".Rdata"), sep = "/")
            

            ### Analalyze tube T ========================================
          } else if (case$TUBE[cs] == "tubeT") {
            
            tryCatch(
              {
            # Read T template content
            # file_contents <- readLines(paste0(
            #   AUTOFLOW_destination, "/R/",
            #   T_function
            # ))

            # condition to select normal or QC analyzing file
            if(grepl("QC", s, ignore.case = T) & grepl("CHEX", s, ignore.case = T) ) {
              
              warning("Loading QC Analyzing template")
              file_contents <- readLines(paste0(AUTOFLOW_destination, "/R/", 
                                                gsub(".Rmd", "", T_function),
                                                "_QC.Rmd"))
            } else {
              file_contents <- readLines(paste0(AUTOFLOW_destination, "/R/", 
                                                T_function))
            }
            
            
            
            
            
            
            file_contents <- file_contents[1:(which(file_contents == "# DEV PART #####################################################################") - 1)]

            fcs.names <- case$FCS_NAME[cs]

            # T_FCS_PATH <- paste(temp_fcs_folder, T_FCS, sep = "/")


            ## a list of content that need to be updated
            list.content_update <- list(
              `afpath`= c("paste0('~/Developer/AUTOFLOW32', '/')", AUTOFLOW_destination),
              `filename` = c("XXX_FILENAME_T_XXX", fcs.names),
              `filepath` = c("XXX_FILEPATH_T_XXX", paste(temp_fcs_folder, fcs.names, sep = "/")),
              `comppath` = c("XXX_COMP_PATH_T_XXX", T_comp_path_local),
              `compfile` = c("XXX_COMP_FILE_T_XXX", T_comp_filename),
              #`pngfolder` = c("XXX_PNG.PATH_XXX", paste(Results_folder, "png", sep = "/")),
              `path.output` = c("/Users/flowuser/Documents/Debug/OUTPUT", path.output),
              `batchfolder` = c("XXX_BATCHFOLDER_XXX", batch),
              `max_analyzing_number` = c("85678935", max_analyzing_number)
            )
            ## Update T script content ===============
            for (l in seq(length(list.content_update))) {
              file_contents <- gsub(
                x = file_contents,
                pattern = list.content_update[[l]][1],
                replacement = list.content_update[[l]][2]
              )
            }

            ### save the updated T analysis file
            temp_t_script_path <- paste(temp_scripts_folder, paste(batch, case$BF[cs], "T.Rmd", sep = "_"), sep = "/")
            cat(file_contents, file = temp_t_script_path, sep = "\n")

            # Run T script ==================
            rmarkdown::render(temp_t_script_path, # output_file = "aa.html",
              output_dir = paste(Results_folder, "html", sep = "/"),  envir = new.env(),
              clean = T
            )
            
            
          }, error = function(e) {
            s <- paste0("FAILED_", s, "TubeT.txt")
            sink(paste0(paste(Results_folder, "reports", sep = "/"), "/", "FAILED_", s, ".txt"))
            cat("ERROR :", conditionMessage(e), "\n")
            sink()
          }
    )

            data_path <- paste(Results_folder, "data", paste0("l.T.data.", StrExtract(fcs.names, ".fcs", 1), ".Rdata"), sep = "/")
            # names(T.data_path) <- fcs.names
          }

          l.case.data[[cs]] <- c(case$TUBE[cs], data_path)
        }

        names(l.case.data) <- case$FCS_NAME

        ### Plotting ====================

        # No fancy theme
        tissue <- case$TISSUE[1]
        # if (tissue == "LN") {
        #   Theme <- "bootstrap"
        # } else if (tissue == "BM") {
        #   Theme <- "yeti"
        # } else if (tissue == "PB") {
        #   Theme <- "united"
        # } else {
        #   Theme <- "default"
        # }
        Theme <- "bootstrap"
        title <- s#paste(case$BF[1], tissue, case$PATIENT[1], case$BCCA[1], sep = "_")


        # lets break down the report into different part

        # report_contents <- readLines(paste0(
        #   AUTOFLOW_destination, "/R/",
        #   REPORT_function
        # )) # readLines(paste(Temp_scripts_Folder, REPORT_function, sep = "/"))
        
        # condition to select normal or QC analyzing file
        if(grepl("QC", s, ignore.case = T) & grepl("CHEX", s, ignore.case = T) ) {
          
          warning("Loading QC Analyzing template")
          report_contents <- readLines(paste0(AUTOFLOW_destination, "/R/", 
                                            gsub(".Rmd", "", REPORT_function),
                                            "_QC.Rmd"))
        } else {
          report_contents <- readLines(paste0(AUTOFLOW_destination, "/R/", 
                                              REPORT_function))
        }
        

        b_starting_line <- which(report_contents == "```{r B Start, include=FALSE}")
        t_starting_line <- which(report_contents == "```{r T Start, include=FALSE}")
        b_end_line <- t_starting_line - 1
        t_end_line <- which(report_contents == "About {.storyboard}") - 1

        head <- report_contents[1:(b_starting_line - 1)]
        footer <- report_contents[(t_end_line + 1):length(report_contents)]
        
        footer <- gsub(
          x = footer,
          pattern = "XXXXXX_STATS_XXXXXXX",
          replacement = paste0(stats_folder, "/", title)
        )

        l.updates <- lapply(seq(length(l.case.data)), function(i) {
          x <- l.case.data[[i]]
          if (is.null(x)) {
            NULL
          } else if (x[1] == "tubeB") {
            
            report_body_b <- report_contents[b_starting_line:b_end_line]
            
            list.content_update <- list(
              `t.datapth` = c("B.DATAPATH", x[2]),
              `t.fcsname` = c("BFCSNAME", case$FCS_NAME[i])
            )
            
            
            #### Update B report script body ===============
            for (l in seq(length(list.content_update))) {
              report_body_b <- gsub(
                x = report_body_b,
                pattern = list.content_update[[l]][1],
                replacement = list.content_update[[l]][2]
              )
            }
            
            report_body_b            
            
            

          } else if (x[1] == "tubeT") { # update T code

            report_body_t <- report_contents[t_starting_line:t_end_line]


            list.content_update <- list(
              `t.datapth` = c("T.DATAPATH", x[2]),
              `t.fcsname` = c("TFCSNAME", case$FCS_NAME[i])
            )


            #### Update T report script body ===============
            for (l in seq(length(list.content_update))) {
              report_body_t <- gsub(
                x = report_body_t,
                pattern = list.content_update[[l]][1],
                replacement = list.content_update[[l]][2]
              )
            }

            report_body_t
          }
        })

        report_contents <- c(head, do.call(c, l.updates), footer)

        ## a list of content that need to be updated
        list.content_update <- list(
          `theme` = c("theme: default", paste("theme:", Theme, sep = " ")),
          `tissue` = c("BCCA Clinical Flow Lab", paste0("Tissue: ", tissue)),
          `title` = c("TITLENAME", title),
          `dpi` = c("DPI", DPI)
        )

        ## Update metadata ===============
        for (l in seq(length(list.content_update))) {
          report_contents <- gsub(
            x = report_contents,
            pattern = list.content_update[[l]][1],
            replacement = list.content_update[[l]][2]
          )
        }

        report_path <- paste0(temp_scripts_folder, "/", title, ".Rmd")
        cat(report_contents, file = report_path, sep = "\n")



        ## Generate Report ===============
        if (is.null(report_only_folder)) {
          rmarkdown::render(report_path,
            output_dir = paste(Results_folder, "reports", sep = "/"), # Output dir make sure the files go to the reports folder
            clean = TRUE
          )
        } else {
          rmarkdown::render(report_path,
            output_dir = report_only_folder, # Output dir make sure the files go to the reports folder
            clean = TRUE
          )
        }

      
        s# <- paste0(s, ".html")
      },
      error = function(e) {
        s <- paste0("FAILED_", s)
        sink(paste0(paste(Results_folder, "reports", sep = "/"), "/", "FAILED_", s, ".txt"))
        cat("ERROR :", conditionMessage(e), "\n")
        sink()
        return(s)
      }
    )
    
    Report_filenames <- c(Report_filenames, log)
  }
  
  
  # Remove files from temp folder
  if(keeptimepfiles == F) {
    file.remove(list.files(paste(Results_folder, "fcs", sep = "/"), full.names = TRUE))
    file.remove(list.files(paste(Results_folder, "html", sep = "/"), full.names = TRUE))
    file.remove(list.files(paste(Results_folder, "png", sep = "/"), full.names = TRUE))
  }

  
  
  
  message(paste0(" Finish analysis for the following samples:\n ", paste(names(file.list), collapse = "\n ")))
  return(Report_filenames)
}







f.quickumap <- function(obj, plotting_clusters, show_BG = F, show_clustername = T) {
  lim <- c(0, 1024)
  umap_cluster <- exprs(obj)[, c("UMAP1", "UMAP2", "Cluster")]
  if (show_clustername == T) {
    l <- split(data.frame(umap_cluster)[, 1:2], umap_cluster[, "Cluster"])
    l <- lapply(l, function(x) {
      apply(x, 2, median)
    })
  }

  if (show_BG == T) {
    dat <- umap_cluster
    dat[, "Cluster"][!umap_cluster[, "Cluster"] %in% plotting_clusters] <- 63
  } else {
    dat <- umap_cluster[umap_cluster[, "Cluster"] %in% plotting_clusters, ]
  }


  rgbwt <- scatter_points_rgbwt(dat[, 1:2],
    out_size = c(256, 256), xlim = lim, ylim = lim,
    palette = col2rgb(hcl_single_color, alpha = T), map = dat[, "Cluster"]
  )
  rstr <- rgba_int_to_raster(rgbwt_to_rgba_int(rgbwt))
  par(pty = "s", bty = "n", mar = c(0, 0, 0, 0), bg = "#F5F5F5") # make the plot square
  # plot(rstr,  interpolate=F)

  plot(c(),
    xlim = lim, xlab = "", xaxt = "n", # remove ticks and label
    ylim = lim, ylab = "", yaxt = "n",
    # main = title, cex.main = 1
  )
  # title(title, adj = 0, line = 0.2)
  rasterImage(rstr,
    xleft = lim[1],
    xright = lim[2],
    ybottom = lim[1],
    ytop = lim[2]
  )



  if (show_clustername == T) {
    lapply(names(l), function(i) {
      shadowtext(l[[i]][1], l[[i]][2], i, cex = 1.1, col = "black", bg = "white", r = 0.2, font = 2)
    })
  }
}








#' Map Numeric Values to Color Gradients
#'
#' This function maps a vector of numeric values to a color gradient using color palettes from the `viridis` package. It supports multiple gradients, including `magma`, `inferno`, `plasma`, `viridis`, and `cividis`. Robust normalization based on percentiles is used to handle outliers.
#'
#' @param numbers A numeric vector containing the values to be mapped to colors.
#' @param palette A character string specifying the color palette to use. Options are `"magma"`, `"inferno"`, `"plasma"`, `"viridis"`, and `"cividis"`. Default is `"magma"`.
#' @param lower_percentile A numeric value between 0 and 1 specifying the lower percentile for robust normalization. Default is `0.01` (1st percentile).
#' @param upper_percentile A numeric value between 0 and 1 specifying the upper percentile for robust normalization. Default is `0.99` (99th percentile).
#' @param breaks An integer specifying the number of color bins in the gradient. Default is `256`, which is the number of colors generated by the `viridis` palettes.
#'
#' @return A character vector of hexadecimal color codes corresponding to the input numbers.
#'
#' @examples
#' # Generate a vector of numbers
#' set.seed(123)
#' numbers <- c(runif(1e6, 1, 100), 1000, -500) # 1 million numbers + 2 outliers
#'
#' # Map numbers to colors using the inferno palette
#' colors <- number_to_color(numbers, palette = "inferno")
#'
#' # Check the first few colors
#' head(colors)
#'
#' @seealso \code{\link[viridis]{magma}}, \code{\link[viridis]{inferno}}, \code{\link[viridis]{plasma}}, \code{\link[viridis]{viridis}}, \code{\link[viridis]{cividis}} for the color palettes used in this function.
#' @export
number_to_color <- function(numbers, palette = "magma", lower_percentile = 0.01, upper_percentile = 0.99, breaks = 256) {
  # Validate the palette input
  if (!palette %in% c("magma", "inferno", "plasma", "viridis", "cividis")) {
    stop("Invalid palette. Choose from 'magma', 'inferno', 'plasma', 'viridis', or 'cividis'.")
  }

  # Compute robust lower and upper bounds using percentiles
  lower_bound <- quantile(numbers, probs = lower_percentile, na.rm = TRUE)
  upper_bound <- quantile(numbers, probs = upper_percentile, na.rm = TRUE)

  # Normalize the numbers to the range [0, 1] using robust bounds
  normalized_numbers <- (numbers - lower_bound) / (upper_bound - lower_bound)

  # Clip values outside the range [0, 1] to avoid extrapolation
  normalized_numbers <- pmin(pmax(normalized_numbers, 0), 1)

  # Map the normalized numbers to the selected color palette
  color_function <- switch(palette,
    magma = viridis::magma,
    inferno = viridis::inferno,
    plasma = viridis::plasma,
    viridis = viridis::viridis,
    cividis = viridis::cividis
  )

  colors <- color_function(breaks)[as.numeric(cut(normalized_numbers, breaks = breaks))]

  return(colors)
}






#' Compute Jaccard coefficient between nearest-neighbor sets.
#'
#' Weights of both i->j and j->i are recorded if they have an intersection. In this case,
#' w(i->j) should be equal to w(j->i). In some cases, i->j has weights while j<-i has no
#' intersection; only w(i->j) is recorded. This is determined in the code by `if (u > 0)`.
#' In this way, the undirected graph is symmetrized by halving the weight
#' in the code `weights[r, 3] <- u / (2.0 * ncol - u) / 2`.
#'
#' @param idx A numeric matrix where each row represents a node and each column
#'   represents its nearest neighbors (indices of other nodes).  Node indices
#'   should be integers.
#'
#' @return A numeric matrix with three columns:
#'   \itemize{
#'     \item Column 1: Index of the first node (i).
#'     \item Column 2: Index of the second node (j).
#'     \item Column 3: Jaccard coefficient between the neighbor sets of i and j.
#'   }
#'   The matrix contains only the pairs (i, j) for which the Jaccard coefficient
#'   is greater than 0.
#'
#' @examples
#' idx_matrix <- matrix(c(1, 2, 3, 2, 1, 4, 3, 4, 1), nrow = 3, byrow = TRUE)
#' result <- jaccard_coeff_r(idx_matrix)
#' print(result)
#'
#' idx_matrix2 <- matrix(c(1, 2, 3, 2, 1, 4, 3, 4, 5), nrow = 3, byrow = TRUE)
#' result2 <- jaccard_coeff_r(idx_matrix2)
#' print(result2)
#'
#' idx_matrix3 <- matrix(c(1, 2, 3, 2, 1, 3, 3, 1, 1), nrow = 3, byrow = TRUE)
#' result3 <- jaccard_coeff_r(idx_matrix3)
#' print(result3)
#'
f.jaccard_coeff <- function(idx) {
  nrow <- nrow(idx)
  ncol <- ncol(idx)
  weights <- matrix(NA, nrow = nrow * ncol, ncol = 3)
  r <- 1

  for (i in 1:nrow) {
    for (j in 1:ncol) {
      k <- idx[i, j]
      nodei <- idx[i, ]
      nodej <- idx[k, ]
      u <- length(intersect(nodei, nodej))

      if (u > 0) {
        weights[r, 1] <- i
        weights[r, 2] <- k
        weights[r, 3] <- u / (2.0 * ncol - u) / 2
        r <- r + 1
      }
    }
  }

  weights <- weights[1:(r - 1), , drop = FALSE]
  return(weights)
}









f.SOM_clustering <- function(obj, clustering.markers,
                             SOM.grid = 6,
                             repeatSOM = 3,
                             k.meta = 6,
                             meta.method = "SOM",
                             show.SOM.plotS = T,
                             show.META.plotS = T) {
  cat(identifier(obj), ":", nrow(obj), "B cells \n")
  channels.ind <- Find.markers(obj, clustering.markers)

  df <- exprs(obj)


  if (nrow(df) <= 20) {
    cat("less than 20 cells, no clustering performed...")
    SOM <- rep(1, nrow(df))
    SOM_META <- data.frame(SOM = SOM, META = SOM)


    if (show.SOM.plotS == T) {
      par(mfrow = c(5, 5))
      flowPlot(df, c("YG585-A", "R675-A"), density.overlay = F, xlim = c(0, 4.5), ylim = c(0, 4.5), title = "<20 cells no clustering")
    }

    if (show.META.plotS == T) {
      par(mfrow = c(5, 5))


      flowPlot(df, c("YG585-A", "R675-A"), density.overlay = F, xlim = c(0, 4.5), ylim = c(0, 4.5), title = "<20 cells no clustering")
    }
  } else {
    # transform FSC SSC
    a <- c <- 0
    b <- 0.004
    df[, c("FSC-A", "SSC-A")] <- apply(df[, c("FSC-A", "SSC-A")], 2, function(x) {
      asinh(a + b * x) + c
    })
    df[, c("FSC-A", "SSC-A")] <- rescale(df[, c("FSC-A", "SSC-A")], c(-0.5, 4))


    # Anonymize the light chain
    # IgL<-fb_results[,colnames(fb_results)%in%c("Eu151Di","Gd160Di")]
    IgL <- df[, c(grep("R675-A", colnames(df)), grep("YG585", colnames(df)))]

    # Set limits
    if (nrow(IgL) > 500) {
      IgL <- apply(IgL, 2, function(x) {
        q <- quantile(x, c(0.005, 0.995))
        x <- pmax(x, q[1])
        x <- pmin(x, q[2])
        # x[x >= q[2]] <- q[2]
        # x[x <= q[1]] <- q[1]
        x
      })
    }

    IgL <- apply(IgL, 2, function(x) {
      x <- rescale(x, c(min(IgL), max(IgL)))
    })


    KL <- apply(IgL, 1, max)

    df <- cbind(df, KL) # remove KL

    # df <- df[, !colnames(df) %in% c("R675-A", "YG585-A")]
    clustering.index <- c(1, 4, channels.ind, ncol(df))

    cat(
      "Clustering with: \n",
      clustering.markers[1:10], "\n",
      clustering.markers[11:20], "\n",
      clustering.markers[21:length(clustering.markers)], "\n",
      "plus FSC, SSC and anonymized light chains...\n"
    )

    Seeds <- c(23L, 67L, 86L, 1L, 27L, 68L, 8L, 92L, 18L, 91L, 133L, 783L, 516L, 2374L, 925L, 7L, 37L, 123L, 49L, 13L)[1:repeatSOM]


    if (nrow(df) <= SOM.grid^2) {
      SOM.grid <- floor(sqrt(10))
    }



    t.som <- system.time(
      l.MAP <- lapply(Seeds, function(s) {
        set.seed(s)
        MAP <- FlowSOM::SOM(df[, clustering.index],
          xdim = SOM.grid, ydim = SOM.grid # ,
          # rlen = 20, mst = 2, alpha = c(0.05, 0.01),#, # radius = stats::quantile(nhbrdist, 0.67) * c(1, 0),
          # init = FALSE, distf = 4, silent = FALSE, codes = NULL#, importance = Markers_importance # weighted markers
        )
        # SOM <-MAP$mapping[,1]
      })
    )
    names(l.MAP) <- Seeds

    cat("DONE ~", t.som[3], "s\n", " Perform", SOM.grid, "X", SOM.grid, "SOM for", length(Seeds), "repeat(s)...")


    l.SOM <- lapply(l.MAP, function(x) {
      x$mapping[, 1]
    })

    if (length(Seeds) > 1) {
      ari.scores <- do.call(rbind, lapply(1:length(names(l.SOM)), function(i) {
        ari.scores <- unlist(lapply(setdiff(1:length(names(l.SOM)), i), function(j) {
          aricode:::ARI(l.SOM[[i]], l.SOM[[j]])
        }))
      }))
      MAP <- l.MAP[[order(apply(ari.scores, 2, median), decreasing = T)[1]]]
    } else {
      MAP <- l.MAP[[1]]
    }

    k.meta <- k.meta
    if (meta.method == "SOM") {
      metaClustering <- FlowSOM::metaClustering_consensus(MAP$codes, k = k.meta, seed = 23)
      dic <- data.frame(SOM = 1:SOM.grid^2, META = metaClustering)
    } else { # h clustering based.

      l.df <- split(as.data.frame(df), f = MAP$mapping[, 1])
      df.median <- do.call(rbind, lapply(l.df, function(x) {
        apply(x[, clustering.index], 2, median)
      }))

      adj.matrix <- 1 - cor(t(df.median))
      hc <- hclust(dist(adj.matrix)) # this is not quite right, but using dist instead of as.dist giving better results
      # plot(hc)


      dic <- data.frame(SOM = 1:SOM.grid^2, META = cutree(hc, k = k.meta))
    }

    SOM <- MAP$mapping[, 1]
    SOM_META <- data.frame(SOM = MAP$mapping[, 1], META = dic$META[match(SOM, dic$SOM)])



    if (show.SOM.plotS == T) {
      l.df <- split(as.data.frame(df), f = SOM_META[, 1])

      par(mfrow = c(5, 5))
      for (i in names(l.df)) {
        flowPlot(l.df[[i]], c("YG585-A", "R675-A"), density.overlay = F, xlim = c(0, 4.5), ylim = c(0, 4.5), title = paste0("SOM ", i))
      }
    }



    ReOrder <- as.numeric(names(sort(table(SOM_META[, 2]), decreasing = T)))

    ReOrder <- data.frame(ReOrder = sort(unique(SOM_META[, 2])), Ori = ReOrder)
    SOM_META[, 2] <- ReOrder$ReOrder[match(SOM_META[, 2], ReOrder$Ori)]


    if (show.META.plotS == T) {
      l.df <- split(as.data.frame(df), f = SOM_META[, 2])
      par(mfrow = c(5, 5))

      for (i in names(l.df)) {
        flowPlot(l.df[[i]], c("YG585-A", "R675-A"), density.overlay = F, xlim = c(0, 4.5), ylim = c(0, 4.5), title = paste0("META ", i))
      }
    }
  }


  # print(apply(SOM_META, 2, table))

  return(SOM_META)
}




#' gate kappa lambda


# x <- cytoframe_to_flowFrame(fs[[1]]); show.plots = T
f.KL <- function(x, show.plots = T) {
  
  dummy_filter <- matrix(rep(-Inf, 6), nrow = 3)
  colnames(dummy_filter) <- c("YG585-A", "R675-A")



  if (nrow(x) <=20) {
    # if (show.plots == T) {
    #   flowPlot(x, channels = c("YG585-A", "R675-A"), xlim = c(0, 4.5), ylim = c(0, 4.5), density.overlay = F)
    #   # lines(IgK_gate@boundaries, col = "blue", lty = 1, lwd = 2)
    #   # lines(IgL_gate@boundaries, col = "green", lty = 1, lwd = 2)
    #   # lines(IgNULL_gate@boundaries, col = "red", lty = 1, lwd = 2)
    #   # lines(IgDP_gate@boundaries, col = "yellow", lty = 1, lwd = 2)
    # }


    IgK_gate <- polygonGate(filterId = "IgK", .gate = dummy_filter)
    IgL_gate <- polygonGate(filterId = "IgL", .gate = dummy_filter)
    IgNULL_gate <- polygonGate(filterId = "IgNULL", .gate = dummy_filter)
    IgDP_gate <- polygonGate(filterId = "IgDP", .gate = dummy_filter)



    return(list = c(
      IgK = IgK_gate,
      IgL = IgL_gate,
      IgNULL = IgNULL_gate,
      IgDP = IgDP_gate
    ))
  } else {
    pks_K <- getPeaks(x, "R675-A")
    pks_L <- getPeaks(x, "YG585-A")


    ct_K <- f.roughGate(x, "R675-A", thrd = 1.5)
    ct_L <- f.roughGate(x, "YG585-A", thrd = 1.5)


    if (show.plots == T) {
      par(mfrow = c(1, 4))

      flowPlot(x, channels = c("YG585-A", "R675-A"), xlim = c(0, 4.5), ylim = c(0, 4.5), title = identifier(x))
      abline(h = ct_K, v = ct_L, col = "black")
    }

    rot.theta <- 0.7 # atan(slope) *1.3  # Angle to rotate by

    # rotate the data
    rot <- rotate.data(x, c("YG585-A", "R675-A"), theta = rot.theta)$data


    #ct_k.rot <- ct_l.rot <- f.roughGate(rot,  "R675-A", thrd = 0)
    ct_k.rot <- f.roughGate(rot,  "R675-A", thrd = 0.5)
    ct_l.rot <- f.roughGate(rot,  "R675-A", thrd = -0.5)
    
    if(ct_l.rot < -1) {ct_l.rot <- ct_k.rot}

    
    

    if (show.plots == T) {
      flowPlot(rot, channels = c("YG585-A", "R675-A"), xlim = c(1, 5), ylim = c(-2, 2))
      abline(h = ct_k.rot, col = "red")
      abline(h = ct_l.rot, col = "green")
    }



    # Perform gating
    cp_K.rot <- flowDensity(rot, c("YG585-A", "R675-A"), position = c(NA, T), gates = c(NA, ct_k.rot))
    cp_L.rot <- flowDensity(rot, c("YG585-A", "R675-A"), position = c(NA, F), gates = c(NA, ct_l.rot))
    
    
    
    
    
    
    
    ## IgNULL
    cp_IgAb.rot <- flowDensity(rot, c("YG585-A", "R675-A"), position = c(NA, F), gates = c(NA, ct_k.rot))
    cp_IgAb.rot <- flowDensity(cp_IgAb.rot, c("YG585-A", "R675-A"), position = c(NA, T), gates = c(NA, ct_l.rot))


    # IgDP
    # pks.null.rot <- getPeaks(cp_IgNULL.rot, "YG585-A", tinypeak.removal = 1/50)
    if (cp_IgAb.rot@cell.count > 5) {
      ct_dp <- f.roughGate(cp_IgAb.rot, "YG585-A", thrd = 3)
    } else {
      ct_dp <- f.roughGate(rot, "YG585-A", thrd = 3)
    }
    cp_IgNULL.rot <- flowDensity(cp_IgAb.rot, c("YG585-A", "R675-A"), position = c(F, NA), gates = c(ct_dp, NA))
    cp_IgDP.rot <- flowDensity(cp_IgAb.rot, c("YG585-A", "R675-A"), position = c(T, NA), gates = c(ct_dp, NA))


    if (show.plots == T) {
      flowPlot(rot, channels = c("YG585-A", "R675-A"), xlim = c(1, 5), ylim = c(-2, 2))
      abline(h = ct_k.rot, col = "red")
      abline(h = ct_l.rot, col = "green")
      abline(v = ct_dp, col = "blue")
    }



    # cp_K <-rotate.fd(cp_K.rot, angle = rot.theta)
    # cp_L<-rotate.fd(cp_L.rot, angle = rot.theta)
    # cp_IgNULL<-rotate.fd(cp_IgNULL.rot, angle = rot.theta)
    # cp_IgDP<-rotate.fd(cp_IgDP.rot, angle = rot.theta)
    #


    if (cp_K.rot@cell.count < 3) {
      IgK_gate <- polygonGate(filterId = "IgK", .gate = dummy_filter)
    } else {
      cp_K <- rotate.fd(cp_K.rot, angle = rot.theta)
      cp_K@filter <- f.GatePrettify(cp_K@filter, ul =T, ur =T, ll =T)
      
      IgK_gate <- polygonGate(filterId = "IgK", .gate = cp_K@filter)
      
    }

    if (cp_L.rot@cell.count < 3) {
      IgL_gate <- polygonGate(filterId = "IgL", .gate = dummy_filter)
    } else {
      cp_L <- rotate.fd(cp_L.rot, angle = rot.theta)
      cp_L@filter <- f.GatePrettify(cp_L@filter, ll =T, ur =T, lr =T)
      
      IgL_gate <- polygonGate(filterId = "IgL", .gate = cp_L@filter)
    }

    if (cp_IgNULL.rot@cell.count < 3) {
      IgNULL_gate <- polygonGate(filterId = "IgNULL", .gate = dummy_filter)
    } else {
      cp_IgNULL <- rotate.fd(cp_IgNULL.rot, angle = rot.theta)
      IgNULL_gate <- polygonGate(filterId = "IgNULL", .gate = cp_IgNULL@filter)
    }

    if (cp_IgDP.rot@cell.count < 3) {
      IgDP_gate <- polygonGate(filterId = "IgDP", .gate = dummy_filter)
    } else {
      cp_IgDP <- rotate.fd(cp_IgDP.rot, angle = rot.theta)
      IgDP_gate <- polygonGate(filterId = "IgDP", .gate = cp_IgDP@filter)
    }






    # cp_K <- f.remove_pop_outliers(cp_K, c("YG585-A", "R675-A"), axis.x = T, axis.y = T)
    # cp_L <- f.remove_pop_outliers(cp_L, c("YG585-A", "R675-A"), axis.x = T, axis.y = T)


    if (show.plots == T) {
      flowPlot(x, channels = c("YG585-A", "R675-A"), xlim = c(0, 4.5), ylim = c(0, 4.5))
      lines(IgK_gate@boundaries, col = "blue", lty = 1, lwd = 2)
      lines(IgL_gate@boundaries, col = "green", lty = 1, lwd = 2)
      lines(IgNULL_gate@boundaries, col = "red", lty = 1, lwd = 2)
      lines(IgDP_gate@boundaries, col = "yellow", lty = 1, lwd = 2)
    }

    return(list = c(
      IgK = IgK_gate,
      IgL = IgL_gate,
      IgNULL = IgNULL_gate,
      IgDP = IgDP_gate
    ))
  }
}



















intersect_polygon_horizontal <- function(polygon, y_val = 0) {
  # Extract x and y coordinates from the polygon matrix
  x <- polygon[, 1]
  y <- polygon[, 2]
  n <- length(x)
  intersections <- list()
  
  for (i in 1:n) {
    # Current and next vertex (with wrap-around for the last edge)
    x1 <- x[i]
    y1 <- y[i]
    x2 <- x[i %% n + 1]
    y2 <- y[i %% n + 1]
    
    # Case 1: Edge is horizontal and lies exactly on y_val
    if (y1 == y2 && y1 == y_val) {
      intersections <- c(intersections, list(c(x1, y_val), c(x2, y_val)))
    } else {
      # Check if the edge crosses y_val
      if ((y1 <= y_val && y_val <= y2) || (y2 <= y_val && y_val <= y1)) {
        # Avoid division by zero (non-horizontal edge guaranteed here)
        t <- (y_val - y1) / (y2 - y1)
        x_intersect <- x1 + t * (x2 - x1)
        intersections <- c(intersections, list(c(x_intersect, y_val)))
      }
    }
  }
  
  # Remove duplicates and handle edge cases
  if (length(intersections) == 0) {
    return(matrix(numeric(0), ncol = 2, dimnames = list(NULL, colnames(polygon))))
  }
  
  unique_intersections <- unique(do.call(rbind, intersections))
  sorted_intersections <- unique_intersections[order(unique_intersections[, 1]), , drop = FALSE]
  colnames(sorted_intersections) <- colnames(polygon)
  sorted_intersections
}

intersect_polygon_line <- function(polygon, y_val = NULL, x_val = NULL) {
  # Validate input: exactly one of y_val or x_val must be specified
  if (sum(c(!is.null(y_val), !is.null(x_val))) != 1) {
    stop("Specify exactly one of: y_val or x_val")
  }
  
  # Extract coordinates
  x <- polygon[, 1]
  y <- polygon[, 2]
  n <- length(x)
  intersections <- list()
  
  if (!is.null(y_val)) {
    # Horizontal line intersection logic
    for (i in 1:n) {
      x1 <- x[i]
      y1 <- y[i]
      x2 <- x[i %% n + 1]
      y2 <- y[i %% n + 1]
      
      # Case 1: Entire edge is exactly on the horizontal line
      if (y1 == y_val && y2 == y_val) {
        intersections <- c(intersections, list(c(x1, y_val), c(x2, y_val)))
        next
      }
      
      # Case 2: Edge crosses the horizontal line
      if ((y1 < y_val && y_val <= y2) || (y2 < y_val && y_val <= y1)) {
        t <- (y_val - y1) / (y2 - y1)
        x_int <- x1 + t * (x2 - x1)
        intersections <- c(intersections, list(c(x_int, y_val)))
      }
    }
  } else {
    # Vertical line intersection logic
    for (i in 1:n) {
      x1 <- x[i]
      y1 <- y[i]
      x2 <- x[i %% n + 1]
      y2 <- y[i %% n + 1]
      
      # Case 1: Entire edge is exactly on the vertical line
      if (x1 == x_val && x2 == x_val) {
        intersections <- c(intersections, list(c(x_val, y1), c(x_val, y2)))
        next
      }
      
      # Case 2: Edge crosses the vertical line
      if ((x1 < x_val && x_val <= x2) || (x2 < x_val && x_val <= x1)) {
        t <- (x_val - x1) / (x2 - x1)
        y_int <- y1 + t * (y2 - y1)
        intersections <- c(intersections, list(c(x_val, y_int)))
      }
    }
  }
  
  # Handle empty case
  if (length(intersections) == 0) {
    return(matrix(numeric(0), 
                  ncol = 2, 
                  dimnames = list(NULL, colnames(polygon))))
  }
  
  # Process results
  unique_points <- unique(do.call(rbind, intersections))
  
  # Sort based on line type
  if (!is.null(y_val)) {
    sorted_points <- unique_points[order(unique_points[, 1]), , drop = FALSE]
  } else {
    sorted_points <- unique_points[order(unique_points[, 2]), , drop = FALSE]
  }
  
  colnames(sorted_points) <- colnames(polygon)
  sorted_points
}








#' Predefined color palettes for flow cytometry populations
#' @export
autoflow_colors <- list(
  `autoflow_30c_T` = c(
    `CD8pCD57p` = "#0000A6",
    `abT/CD4pCD8n` = "#FF4444",
    `abT/CD4nCD8p` = "#FF4444",
    #`CD3pCD2p` = "#BC23FF",
    `CD56pCD3n` = "#BE26F6", #"#3A2465", 
    `B` = "#0776f5", # "#08306B",
    `CD3pCD19p` = "#FFDBE5",
    `Lymph` = "#8AD466",
    `Mono` = "#FFC24F",
    `Live` = "#525252",
    `Others` = "#A75725",
    `root` = "#61615A"
  ),
  `autoflow_30c_B` = c(
    `B_Pop_4` = "#220145",#430a80", 
    `B_Pop_3` = "#0000A6", 
    `B_Pop_2` = "#0776f5", 
    `B_Pop_1` = "#027d9c",
    `B` = "#08306B",
    `T` = "#FF4444",
    `NonTB` = "#8AD466", 
    #`CD3pCD19p` = "#FFDBE5",
    `Lymph` = "#8AD466",
    `Mono` = "#FFC24F",
    `Live` = "#525252",
    `Others` = "#A75725",
    `root` = "#61615A"
  )
)



#' Predefined color palettes for flow cytometry populations
#' @export
autoflow_colors_CLT <- list(
  `autoflow_30c_T` = c(
    `gdT` = "#220145",
    `CD8pCD57p` = "#0000A6",
    `abT/CD4pCD8n` = "#FF4444",
    `abT/CD4nCD8p` = "#FF4444",
    `CD3pCD2p` = "#FF4444",
    `NK` = "#BE26F6", #"#3A2465", 
    `NKT` = "#220145",
    `B` = "#0776f5", # "#08306B",
    `Lymph` = "#8AD466",
    `Mono` = "#FFC24F",
    `Live` = "#525252",
    `Off_scale` = "#A75725",
    `root` = "#61615A"
  ),
  `autoflow_30c_B` = c(
    `B_Pop4` = "#220145",#430a80", 
    `B_Pop3` = "#0000A6", 
    `B_Pop2` = "#0072F5", 
    `B_Pop1` = "#008C7D", #33A6C9
    
    #`Lambda_CD5nCD10n` = "#0000A6", 
    #`Poly_PseudoCD5p` = "#0776f5", 
    #`Poly_Regular` = "#027d9c",
    
    `CLUSTER_B` = "#0776f5",`B` = "#0776f5",#08306B",
    `T` = "#FF4444",
    `CLUSTER_T` = "#FF4444",
    
    `Diag` = "pink",
    `DP` = "pink",
    #`CLUSTER_Baso` = "brown", 
    `Baso` = "black", 
    #`CLUSTER_PDC` = "darkgreen",
    `PDC` = "brown", 
    `NonTB` = "#8AD466", 
    #`CD3pCD19p` = "#FFDBE5",

    `Lymph` = "#8AD466",
    `Mono` = "#FFC24F",

    `Live` = "#525252",
    `Off_scale` = "#A75725",
    `root` = "#61615A"
  ),
  `autoflow_30c_B_QC` = c(

    `B` = "#0776f5",
    `T` = "#FF4444",
    `NonTB` = "darkgreen", 
    `Baso` = "black", 
    #`CLUSTER_PDC` = "darkgreen",
    `PDC` = "brown", 
    `NonTB` = "#8AD466", 
    #`CD3pCD19p` = "#FFDBE5",
    
    `Lymph` = "#8AD466",
    `Mono` = "#FFC24F",
    `Grans` = "#A75725",
    `Live` = "#525252",
    `root` = "#61615A"
  )
)





# classifer function


f.tabpfn <- function(dat, clf) {
  dat <- dat[, colnames(dat) %in% clf$feature_names_in_]
  pred_proba <- clf$predict_proba(dat)
  # colnames(pred_proba) <- clf_Lymph$classes_
  index <- apply(pred_proba, 1, which.max)
  results <- data.frame(cluster = rownames(dat), pred_class = clf$classes_[index])
  results$score <- unlist(lapply(seq(nrow(pred_proba)), function(i) {
    round(pred_proba[i, index[i]], 4)
  }))
  results
}




lighten <- function(color, factor = 0.5) {
  if ((factor > 1) | (factor < 0)) stop("factor needs to be within [0,1]")
  col <- col2rgb(color)
  col <- col + (255 - col)*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}



darken <- function(color, factor = 0.5) {
  if ((factor > 1) | (factor < 0)) stop("factor needs to be within [0,1]")
  nms<- names(color)
  col <- col2rgb(color)
  col <- col * (1 - factor)
  col <- rgb(t(col), maxColorValue = 255)
  names(col) <- nms
  col
}


saturate <- function(color, factor = 0.5) {
  if ((factor > 1) | (factor < 0)) stop("factor needs to be within [0,1]")
  col <- col2rgb(color)
  hsv_col <- rgb2hsv(col)
  hsv_col[2,] <- pmin(hsv_col[2,] * (1 + factor), 1)
  col <- hsv(hsv_col[1,], hsv_col[2,], hsv_col[3,])
  col
}


adjust_saturation <- function(color, factor = 0) {
  if ((factor > 1) | (factor < -1)) stop("factor needs to be within [-1,1]")
  col <- col2rgb(color)
  hsv_col <- rgb2hsv(col)
  saturation <- hsv_col[2,]
  if (factor >= 0) {
    hsv_col[2,] <- pmin(saturation * (1 + factor), 1) # Increase saturation
  } else {
    hsv_col[2,] <- pmax(saturation * (1 + factor), 0) # Decrease saturation
  }
  col <- hsv(hsv_col[1,], hsv_col[2,], hsv_col[3,])
  names(col) <- names (color)
  col
}


PlotFlowFrame <- function (obj, channels, colormap_marker = NA, pop_idex.list = NA, 
          plotting_background = FALSE, title, title.font.size, subtitle, 
          subtitle.color, comp = FALSE, spillovermatix = NA, min.xy, 
          max.xy, lowerlimit.quantile = 3e-04, showticks = TRUE, shownumbers = FALSE, 
          showlabels = TRUE, plot.margins = c(2, 2, 2, 2), res = NA, 
          dens_res, population_colors = Color103, auto_assign_color = TRUE, 
          background_color = "#525252", margin_color, rgba_blending = TRUE, 
          interpolate = FALSE, rare_on_top = TRUE, highlight_pops = NA, 
          highlight.alpha = 0.8, transformation = "logicle", tr.parameters = NA, 
          print_estimated_parameters = FALSE, default_b = 400, log_fsc = FALSE, 
          log_ssc = FALSE, label_size = 0.8, color_pops = TRUE, gates, 
          gates_color = "grey", old_method = TRUE, density.overlay = FALSE, 
          max_points = 5e+05, fulllength_label = T) 
{
  if (class(obj) == "flowFrame") {
  }
  else if (class(obj) == "cytoframe") {
    obj <- cytoframe_to_flowFrame(obj)
  }
  else {
    message("GatingSet or flowSet, only plot the 1st element")
    if (class(obj) == "GatingSet") {
      obj <- gs_pop_get_data(obj, "root")
    }
    if (class(obj) == "flowSet") {
      obj <- obj
    }
    obj <- obj[[1]]
  }
  if (missing(title)) {
    title = identifier(obj)
  }
  channels <- f.parseChannels(obj, channels)
  if (length(channels) < 2) {
    stop("At least one marker is not presented in the data!")
  }
  p_chnames <- parameters(obj)[["name"]][channels]
  p_markers <- as.vector(parameters(obj)[["desc"]])[channels]
  if (length(p_chnames) == 0 || is.null(p_chnames)) {
    stop("Channel names could not be determined. Check that the channels parameter is valid.")
  }
  if (length(p_markers) == 0) {
    p_markers <- rep("", length(p_chnames))
  }
  if (sum(grepl("SC|ime", p_markers)) == 2) {
    comp <- FALSE
  }
  p_markers[is.na(p_markers)] <- ""
  if (comp == TRUE) {
    if (sum(grepl("Comp", p_chnames)) == 0) {
      if (is.na(spillovermatix)) {
        spillovermatix <- spillover(obj)
        spillovermatix <- spillovermatix[!unlist(lapply(spillover(obj), 
                                                        is.null))][[1]]
      }
      obj <- compensate(obj, spillovermatix)
      label_prefix <- "Comp-"
    }
    else {
      p_chnames <- gsub("Comp-", "", p_chnames)
      label_prefix <- "Comp-"
    }
  }
  else {
    label_prefix <- ""
  }
  if (missing(pop_idex.list) | length(na.omit(pop_idex.list)) == 
      0) {
  }
  else {
    population_colors_updated <- lapply(seq(length(pop_idex.list)), 
                                        function(x) {
                                          if (!names(pop_idex.list)[x] %in% names(population_colors)) {
                                            if (auto_assign_color == TRUE) {
                                              warning(paste0("No color assined for ", names(pop_idex.list)[x], 
                                                             ", using  color ", Color103[x]))
                                              Color103[x]
                                            }
                                            else {
                                              background_color
                                            }
                                          }
                                          else {
                                            population_colors[[names(pop_idex.list)[x]]]
                                          }
                                        })
    names(population_colors_updated) <- names(pop_idex.list)
    population_colors_updated <- unlist(population_colors_updated)
    if (plotting_background == TRUE) {
      r <- setdiff(seq(nrow(obj)), unlist(pop_idex.list))
      if (length(r) > 0) {
        pop_idex.list <- c(pop_idex.list, list(root = r))
        population_colors_updated <- c(population_colors_updated, 
                                       background_color)
        names(population_colors_updated) <- names(pop_idex.list)
      }
    }
  }
  dat <- exprs(obj)[, channels]
  if (nrow(dat) > max_points) {
    set.seed(42)
    idx <- sample(nrow(dat), max_points)
    dat <- dat[idx, ]
    if (!missing(pop_idex.list) && length(na.omit(pop_idex.list)) > 
        0) {
      pop_idex.list <- lapply(pop_idex.list, function(x) {
        intersect(x, idx)
      })
    }
  }
  if (!is.null(obj@parameters@data$trans)) {
    trans_parameters_used <- data.frame(row.names = obj@parameters@data$name, 
                                        trans = obj@parameters@data$trans)
  }
  else {
    trans_parameters_used <- NULL
  }
  ranges <- list(c(), c())
  major_tick_loc <- list(c(), c())
  minor_tick_loc <- list(c(), c())
  if (missing(min.xy)) {
    min.xy <- c(NA, NA)
  }
  else {
    if (length(min.xy) == 1) {
      min.xy <- c(min.xy[1], min.xy[1])
    }
  }
  if (missing(max.xy)) {
    max.xy <- c(NA, NA)
  }
  else {
    if (length(max.xy) == 1) {
      max.xy <- c(max.xy[1], max.xy[1])
    }
  }
  for (i in 1:2) {
    if (i > length(p_chnames)) {
      stop(paste("Channel index", i, "is out of bounds. Only", 
                 length(p_chnames), "channels available."))
    }
    if (grepl("SC", p_chnames[i])) {
      if (grepl("SC", p_chnames[i])) {
        min.xy[i] <- 30
      }
      else {
        min.xy[i] <- as.numeric(parameters(obj)[["minRange"]][channels[i]])
      }
      if (is.na(max.xy[i])) {
        plot_ceiling <- as.numeric(parameters(obj)[["maxRange"]][channels[i]])
      }
      else {
        plot_ceiling <- max.xy[i]
      }
      if (grepl("FSC", p_chnames[i])) {
        if (log_fsc == TRUE) {
          mtd <- "log"
        }
        else {
          mtd <- "linear"
        }
      }
      if (grepl("SSC", p_chnames[i])) {
        if (log_ssc == TRUE) {
          min.xy[i] <- 3000
          plot_ceiling <- 262144
          mtd <- "log"
        }
        else {
          mtd <- "linear"
        }
      }
      lim <- c(min.xy[i], plot_ceiling)
      dat[, i][dat[, i] < min.xy[i]] <- min.xy[i]
      dat[, i][dat[, i] > plot_ceiling] <- plot_ceiling
      ranges[[i]] <- f.transFscSscData(lim, method = mtd)
      dat[, i] <- f.transFscSscData(dat[, i], method = mtd)
      if (mtd == "linear") {
        major_tick_loc[[i]] <- major_tick_linear
        minor_tick_loc[[i]] <- minor_tick_linear
      }
      else {
        major_tick_loc[[i]] <- f.transFscSscData(major_tick, 
                                                 method = mtd)
        minor_tick_loc[[i]] <- f.transFscSscData(minor_tick, 
                                                 method = mtd)
      }
    }
    else if (grepl("ime", p_chnames[i])) {
      ranges[[i]] <- c(0, max(dat[, i]))
      time_ticks <- setNames(seq(0, max(dat[, i]), 1000), 
                             seq(0, max(dat[, i]), 1000))
      major_tick_loc[[i]] <- time_ticks
      minor_tick_loc[[i]] <- NA
    }
    else if (grepl("PC|UMAP", p_chnames[i])) {
      if (is.na(min.xy[i])) {
        min.xy[i] <- as.numeric(parameters(obj)[["minRange"]][channels[i]])
      }
      if (is.na(max.xy[i])) {
        plot_ceiling <- as.numeric(parameters(obj)[["maxRange"]][channels[i]])
      }
      else {
        plot_ceiling <- max.xy[i]
      }
      lim <- c(min.xy[i], plot_ceiling)
      dat[, i] <- pmax(dat[, i], lim[1])
      dat[, i] <- pmin(dat[, i], lim[2])
      ranges[[i]] <- f.transFscSscData(lim, method = "linear")
      dat[, i] <- f.transFscSscData(dat[, i], method = "linear")
      major_tick_loc[[i]] <- seq(-10, 10, 1)
      minor_tick_loc[[i]] <- seq(-10, 10, 0.5)
    }
    else {
      if (transformation == "log") {
        min.xy[i] <- 30
      }
      if (is.na(min.xy[i])) {
        plot_floor <- parameters(obj)[["minRange"]][channels[i]]
      }
      else {
        plot_floor <- min.xy[i]
      }
      if (is.na(max.xy[i])) {
        plot_ceiling <- parameters(obj)[["maxRange"]][channels[i]]
      }
      else {
        plot_ceiling <- max.xy[i]
      }
      if (plot_ceiling < 10) {
        if (is.na(min.xy[i]) | (min.xy[i] < -10)) {
          bel <- round(sum(dat[, i] <= plot_floor)/nrow(dat) * 
                         100, 1)
          if (bel >= 10) {
            warning(paste0(bel, " % of events on axis"))
          }
          lim <- c(plot_floor, plot_ceiling)
          dat[, i][dat[, i] < lim[1]] <- lim[1]
        }
        else {
          lim <- c(min.xy[i], plot_ceiling)
          dat[, i][dat[, i] < lim[1]] <- lim[1]
        }
        ranges[[i]] <- lim
        if (is.null(trans_parameters_used)) {
          major_tick_loc[[i]] <- unique(c(ceiling(lim[1]), 
                                          c(0, 1, 2, 3, 4)))
          names(major_tick_loc[[i]]) <- major_tick_loc[[i]]
          minor_tick_loc[[i]] <- major_tick_loc[[i]] + 
            0.5
          names(minor_tick_loc[[i]]) <- minor_tick_loc[[i]]
        }
        else {
          tr.parameter_used <- setNames(trans_parameters_used[p_chnames[i], 
          ], "w")
          major_tick_loc[[i]] <- f.transVectData(major_tick, 
                                                 method = "logicle", tr.parameter_used)
          names(major_tick_loc[[i]]) <- major_tick
          minor_tick_loc[[i]] <- f.transVectData(minor_tick, 
                                                 method = "logicle", tr.parameter_used)
          names(minor_tick_loc[[i]]) <- minor_tick
        }
      }
      else {
        if (length(na.omit(tr.parameters)) == 0) {
          if (transformation == "arcsinh") {
            tr.parameter <- setNames(400, "b")
            message(paste0("No transformation parameters, perform arcsinh transformation with b=400 for ", 
                           p_chnames[i]))
          }
          else if (transformation == "logicle") {
            lgcl <- estimateLogicle(obj, channels = p_chnames[i])
            lgcl_unlisted <- unlist(summary(lgcl))
            if (class(lgcl_unlisted) == "list") {
              tr.parameter <- round(c(lgcl_unlisted[[1]], 
                                      lgcl_unlisted[[3]], lgcl_unlisted[[4]], 
                                      lgcl_unlisted[[6]]), 2)
            }
            else {
              tr.parameter <- round(lgcl_unlisted, 2)
            }
            names(tr.parameter) <- c("a", "m", "t", "w")
            message(paste0("No transformation parameters, perform logicle transformation with estimated w=", 
                           tr.parameter["w"], " for ", p_chnames[i]))
          }
          else {
          }
        }
        else {
          if (length(tr.parameters) == 1) {
            tr.parameters <- c(tr.parameters, tr.parameters)
          }
          if (class(tr.parameters) == "list") {
            tr.parameter <- tr.parameters[[p_chnames[i]]]
          }
          else if (class(tr.parameters) == "numeric") {
            if (transformation == "arcsinh") {
              if (tr.parameters[i] < 100) {
                message("b is too small, set b = 400")
                tr.parameter <- setNames(400, "b")
              }
              else {
                tr.parameter <- setNames(tr.parameters[i], 
                                         "b")
              }
            }
            else if (transformation == "logicle") {
              if (tr.parameters[i] > 5) {
                message("w is too big, set w=0.8")
                tr.parameter <- setNames(0.8, "w")
              }
              else {
                tr.parameter <- setNames(tr.parameters[i], 
                                         "w")
              }
            }
          }
        }
        if (is.na(min.xy[i])) {
          if (is.na(lowerlimit.quantile)) {
            plot_floor <- min(dat[, i])
          }
          else {
            plot_floor <- quantile(dat[, i], lowerlimit.quantile)
          }
          lim <- c(plot_floor, plot_ceiling)
          dat[, i] <- pmax(dat[, i], lim[1])
        }
        else {
          lim <- c(min.xy[i], plot_ceiling)
          dat[, i] <- pmax(dat[, i], lim[1])
        }
        if (print_estimated_parameters) {
          print(unlist(tr.parameter))
        }
        ranges[[i]] <- f.transVectData(lim, method = transformation, 
                                       tr.parameter)
        dat[, i] <- f.transVectData(dat[, i], method = transformation, 
                                    tr.parameter)
        major_tick_loc[[i]] <- f.transVectData(major_tick, 
                                               method = transformation, tr.parameter)
        names(major_tick_loc[[i]]) <- major_tick
        minor_tick_loc[[i]] <- f.transVectData(minor_tick, 
                                               method = transformation, tr.parameter)
        names(minor_tick_loc[[i]]) <- minor_tick
      }
    }
  }
  names(ranges) <- p_chnames
  if (is.na(res)) {
    if (missing(pop_idex.list)) {
      cn <- nrow(dat)
    }
    else {
      cn <- length(unlist(pop_idex.list))
    }
    cat("Number of cells to plot:", cn, "\n")
    if (cn >= 50000) {
      res <- 128
    }
    else if (cn < 50000 & cn >= 5000) {
      res <- 96
    }
    else if (cn < 5000 & cn >= 2500) {
      res <- 64
    }
    else if (cn < 2500 & cn >= 200) {
      res <- 64
    }
    else if (cn < 200 & cn >= 50) {
      res <- 64
    }
    else {
      res <- 32
    }
  }
  else {
    cn <- nrow(dat)
    res <- res
  }
  #cat("auto set resolution to:", res, "\n")
  if (missing(dens_res)) {
    dens_res <- res
  }
  if (is.null(color_pops)) {
    color_pops <- FALSE
  }
  if (cn != 0) {
    if (!is.na(colormap_marker)) {
      colormap_idx <- f.parseChannels(obj, colormap_marker)
      color_name <- as.vector(parameters(obj)[["desc"]])[colormap_idx]
      na <- rowSums(is.na(dat)) > 0
      dat <- dat[!na, ]
      mk <- exprs(obj)[!na, colormap_idx]
      col <- number_to_color(mk, palette = "inferno")
      col <- col2rgb(col, alpha = TRUE)
      rgbwt <- scatter_points_rgbwt(dat, out_size = c(res, 
                                                      res), RGBA = col, xlim = ranges[[1]], ylim = ranges[[2]])
      rgbwt <- rgbwt_to_rgba_int(rgbwt)
    }
    else {
      if (length(na.omit(pop_idex.list)) == 0 | color_pops == 
          FALSE) {
        if (length(na.omit(pop_idex.list)) != 0) {
          l.dat <- lapply(pop_idex.list, function(x) {
            dat[x, ]
          })
          dat <- do.call(rbind, l.dat)
        }
        na <- rowSums(is.na(dat)) > 0
        dat <- dat[!na, ]
        if (old_method) {
          rgbwt <- scatter_points_rgbwt(dat, out_size = c(res, 
                                                          res), RGBA = col2rgb("#25252505", alpha = TRUE), 
                                        xlim = ranges[[1]], ylim = ranges[[2]])
          rgbwt[, , 5] <- log1p(1 - rgbwt[, , 5])/log(2) * 
            0.95
          rgbwt[, , 5][rgbwt[, , 5] == 0] <- 1
        }
        else {
          rgbwt <- scatter_points_rgbwt(dat, out_size = c(res, 
                                                          res), xlim = ranges[[1]], ylim = ranges[[2]])
          rgbwt <- scale_toward_white(rgbwt, base_color = c(0, 
                                                            0, 0), whiteness = 0.5)
        }
        rgbwt <- rgbwt_to_rgba_int(rgbwt)
      }
      else {
        l.dat <- lapply(pop_idex.list, function(x) {
          dat[x, , drop = FALSE]
        })
        hl_pop <- names(l.dat)[names(l.dat) %in% 
                                 highlight_pops]
        pop_sizes <- vapply(l.dat, nrow, numeric(1))
        if (rare_on_top == TRUE) {
          l.dat <- l.dat[names(sort(pop_sizes, decreasing = FALSE))]
        }
        if (old_method) {
          l.rgbwt <- lapply(names(l.dat), function(x) {
            dat_subset <- l.dat[[x]]
            color <- population_colors_updated[x]
            color <- paste0(color, "10")
            rgbwt <- scatter_points_rgbwt(dat_subset, 
                                          out_size = c(res, res), RGBA = col2rgb(color, 
                                                                                 alpha = TRUE), xlim = ranges[[1]], ylim = ranges[[2]])
            if (x %in% highlight_pops) {
              rgbwt <- apply_kernel_rgbwt(rgbwt, "square", 
                                          radius = 1)
            }
            rgbwt[, , 5] <- log1p(1 - rgbwt[, , 5])/0.8
            rgbwt[, , 5][rgbwt[, , 5] == 0] <- 1
            if (rgba_blending == TRUE) {
              rgbwt <- rgbwt_to_rgba_float(rgbwt)
            }
            return(rgbwt)
          })
        }
        else {
          l.rgbwt <- lapply(names(l.dat), function(x) {
            

            if (length(hl_pop) > 0) {
              res <- res*2
            } else {}
            
            dat_subset <- l.dat[[x]]
            color <- population_colors_updated[x]
            rgbwt <- scatter_points_rgbwt(dat_subset, 
                                          out_size = c(res, res), xlim = ranges[[1]], 
                                          ylim = ranges[[2]])
            if (x %in% highlight_pops) {
              cn <- nrow(dat_subset)
              #sr_factor <- 10
              if (cn >= 5000) {
                rds <- 1
                
              }
              else if (cn < 5000 & cn >= 2500) {
                rds <- 1
              }
              else if (cn < 2500 & cn >= 200) {
                rds <- 1#2/sr_factor
              }
              else {
                rds <- 1#3/sr_factor
              }
              rgbwt <- apply_kernel_rgbwt(rgbwt, "square", 
                                          radius = rds)
            }
            rgbwt <- scale_toward_white(rgbwt, base_color = c(col2rgb(color, 
                                                                      alpha = FALSE)), whiteness = 0.95)
            if (rgba_blending == TRUE) {
              rgbwt <- rgbwt_to_rgba_float(rgbwt)
            }
            return(rgbwt)
          })
        }
        if (rgba_blending == TRUE) {
          if (old_method) {
            rgbwt <- blend_rgba_float(l.rgbwt)
          }
          else {

            if (length(hl_pop) > 0) {
              hl_index <- which(names(l.dat) %in% hl_pop)
              l.rgbwt <- c(l.rgbwt[hl_index], l.rgbwt[-hl_index])
              rgbwt <- blend_rgba_float3(l.rgbwt, alpha = highlight.alpha)
            }
            else {
              rgbwt <- blend_rgba_float3(l.rgbwt, alpha = 0.6)
            }
          }
          rgbwt <- rgba_float_to_rgba_int(rgbwt)
        }
        else {
          rgbwt <- merge_rgbwt(l.rgbwt)
          rgbwt <- rgbwt_to_rgba_int(rgbwt)
        }
      }
    }
    rstr <- rgba_int_to_raster(rgbwt)
  }
  usr <- c(ranges[[1]], ranges[[2]])
  par(pty = "s", bty = "o", mgp = c(3, 0.5, 0), mar = plot.margins)
  if (!missing(margin_color)) {
    par(bg = margin_color)
  }
  plot(c(), xlim = ranges[[1]], xlab = "", xaxt = "n", ylim = ranges[[2]], 
       ylab = "", yaxt = "n", )
  if (missing(title.font.size)) {
    if (nchar(title) > 30) {
      title.font.size = 0.6
    }
    else {
      title.font.size = 1
    }
  }
  title(title, adj = 0, line = 0.3, cex.main = title.font.size)
  if (missing(subtitle)) {
    if (!is.na(colormap_marker)) {
      title(color_name, adj = 1, line = 0.3, cex.main = title.font.size)
    }
  }
  else {
    if (missing(subtitle.color)) {
      title(subtitle, adj = 1, line = 0.3, cex.main = title.font.size)
    }
    else {
      title(subtitle, adj = 1, line = 0.3, cex.main = title.font.size, 
            col.main = subtitle.color)
    }
  }
  if (!missing(margin_color)) {
    u <- par("usr")
    rect(u[1], u[3], u[2], u[4], col = "white", border = "black")
    par(new = TRUE)
  }
  if (cn != 0) {
    rasterImage(rstr, xleft = usr[1], xright = usr[2], ybottom = usr[3], 
                ytop = usr[4], interpolate = interpolate)
  }
  if (showticks == TRUE) {
    f.axis_ticks(p_chnames, major_tick_loc, minor_tick_loc, 
                 ranges, shownumbers)
  }
  if (showlabels == TRUE) {
    if (showticks == TRUE & shownumbers == TRUE) {
      label_padding <- 1.3
    }
    else if (showticks == TRUE & shownumbers == FALSE) {
      label_padding <- 0.9
    }
    else {
      label_padding <- 0.5
    }
    axis_labels <- p_chnames
    if (fulllength_label == T) {
      axis_labels[!grepl("SC|ime", p_chnames)] <- paste0(label_prefix, 
                                                         p_chnames[!grepl("SC|ime", p_chnames)], "::", 
                                                         p_markers[!grepl("SC|ime", p_chnames)])
    }
    else {
      axis_labels[!grepl("SC|ime", p_chnames)] <- paste0(label_prefix, 
                                                         p_markers[!grepl("SC|ime", p_chnames)])
    }
    if (log_ssc == T) {
      axis_labels[grepl("SSC", p_chnames)] <- paste0(p_chnames[grepl("SSC", 
                                                                     p_chnames)], " LOG")
    }
    if (log_fsc == T) {
      axis_labels[grepl("FSC", p_chnames)] <- paste0(p_chnames[grepl("FSC", 
                                                                     p_chnames)], " LOG")
    }
    f.axis_label(axis_labels, padding = label_padding, font.size = label_size)
  }
  if (!missing(gates)) {
    f.plotGate(gates, xlim = ranges[[1]], ylim = ranges[[2]], 
               log_ssc = log_ssc, col = gates_color)
  }
  if (density.overlay & cn > 10) {
    dens.color <- "red"
    dens.alpha = 0.5
    adjust.dens = 1
    x.dens <- density(dat[, 1], adjust = adjust.dens)
    y.dens <- density(dat[, 2], adjust = adjust.dens)
    x.axis <- x.dens$x
    y.axis <- y.dens$x
    x.dens <- rescale(x.dens$y, ranges[[2]])
    y.dens <- rescale(y.dens$y, ranges[[1]])
    lines(x.axis, x.dens, col = alpha(dens.color, dens.alpha))
    lines(y.dens, y.axis, col = alpha(dens.color, dens.alpha))
  }
}

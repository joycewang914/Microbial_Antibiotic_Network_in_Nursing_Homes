# Import data stored in csv files
mat = read.csv("example_data_Fig1.csv", header = T)

# Check the dimension of the matrix: 
dim(mat)
# 234 rows = 234 patients
# 6 columns = binary colonization data of 6 prevalent resistant gram-negative bacteria (across all visits)
# Patients ordered by number of bacteria they have 

# First character of patient name indicates which facility the patient was from (A - L)

# Assign each facility with a unique colour
cols = c("violet",  "slateblue1", 
         "turquoise2", "skyblue", 
         "orange", "tomato", "violetred", 
         "springgreen2", "palegreen4",
         "wheat2",  "tan3", "brown")

facil_col = structure(cols, names = LETTERS[1:12])
pt_col = sapply(rownames(mat), FUN = function(x){facil_col[substr(x, 1, 1)]})

# heatmap
heatmap(mat, 
        scale = "none", 
        keep.dendro = FALSE, 
        Rowv = NA,
        Colv = NA,
        labRow = "", 
        labCol = colnames(mat),
        ylab = "Patient ID", 
        main = "Organism by patient", 
        breaks = c(-2, -1, 0, 1), 
        col = c('black', 'white', 'lightblue'),
        RowSideColors = pt_col, 
        margins = c(4,6),
        cexRow = 0.1, 
        cexCol = 1.75)

legend("topright", legend = names(facil_col), fill = facil_col, title = "Facility", bty = "n", border = FALSE)

# To visualize the number of colonization by patient, use barplot function (max = 6)
barplot(rowSums(mat), xaxt = "n")


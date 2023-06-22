library(data.table)

# Migration ----------------------------------------------------------------------------------------
params <- HivEstInfTime::GetMigrantParams()

# Recon data set
reconAIDS <- data.table::setDT(haven::read_dta('D:/VirtualBox_Shared/Migrant_test/baseAIDS.dta'))
reconCD4VL <- data.table::setDT(haven::read_dta('D:/VirtualBox_Shared/Migrant_test/baseCD4VL.dta'))

# Create inputs for testing
baseAIDS <- reconAIDS[, 1:18]
isLabelled <- sapply(baseAIDS, haven::is.labelled)
colNames <- names(isLabelled[isLabelled])
baseAIDS[, (colNames) := lapply(.SD, haven::as_factor), .SDcols = colNames]
data.table::setnames(
  baseAIDS,
  c(
    'RecordId', 'Gender', 'Transmission', 'Age', 'MigrantRegionOfOrigin', 'Calendar', 'Art',
    'DateOfArt', 'DateOfHIVDiagnosis', 'DateOfAIDSDiagnosis', 'DateOfArrival', 'DateOfBirth',
    'AtRiskDate', 'U', 'Mig', 'KnownPrePost', 'UniqueId', 'DTime'
  )
)
currentLevels <- levels(baseAIDS$Gender)
newLevels <- c('M', 'F')
levels(baseAIDS$Gender) <- newLevels[match(currentLevels, c('Male', 'Female'))]
baseAIDS[, Ord := data.table::rowid(UniqueId)]

baseCD4VL <- reconCD4VL[, 1:27]
isLabelled <- sapply(baseCD4VL, haven::is.labelled)
colNames <- names(isLabelled[isLabelled])
baseCD4VL[, (colNames) := lapply(.SD, haven::as_factor), .SDcols = colNames]
data.table::setnames(
  baseCD4VL,
  c(
    'RecordId', 'DateOfExam', 'YVar', 'Indi', 'Gender', 'Transmission', 'Age',
    'MigrantRegionOfOrigin', 'Calendar', 'Art', 'DateOfArt', 'DateOfHIVDiagnosis',
    'DateOfAIDSDiagnosis', 'DateOfArrival', 'DateOfBirth', 'AtRiskDate', 'U', 'Mig', 'KnownPrePost',
    'DTime', 'UniqueId', 'Consc', 'Consr', 'CobsTime', 'RobsTime', 'RLogObsTime2', 'Only'
  )
)
currentLevels <- levels(baseCD4VL$Gender)
newLevels <- c('M', 'F')
levels(baseCD4VL$Gender) <- newLevels[match(currentLevels, c('Male', 'Female'))]
baseCD4VL[, Ord := data.table::rowid(UniqueId)]
baseCD4VL[, UniqueId := UniqueId + max(baseAIDS$UniqueId)]
input <- list(
  Data = list(
    AIDS = baseAIDS,
    CD4VL = baseCD4VL
  )
)

# Create test dataset
test <- HivEstInfTime::PredictInf(input, params)
test <- HivEstInfTime::GetStats(predictions = test, input)

# Reconcile
recon <- rbind(
  unique(reconAIDS[, .(UniqueId = id, ProbPre)]),
  unique(reconCD4VL[, .(UniqueId = id + max(reconAIDS$id), ProbPre)])
)

compare <- merge(
  recon,
  test,
  by = c('UniqueId'),
  suffix = c('.Recon', '.Test'),
  all = TRUE
)
compare[, Diff := ProbPre.Recon - ProbPre.Test]

# Show differences
compare[abs(Diff) > 1e-3]
compare[is.na(ProbPre.Test)]



baseCD4VL2 <- reconCD4VL[, 1:27]
isLabelled <- sapply(baseCD4VL2, haven::is.labelled)
colNames <- names(isLabelled[isLabelled])
baseCD4VL2[, (colNames) := lapply(.SD, haven::as_factor), .SDcols = colNames]
data.table::setnames(
  baseCD4VL2,
  c(
    'RecordId', 'DateOfExam', 'YVar', 'Indi', 'Gender', 'Transmission', 'Age',
    'MigrantRegionOfOrigin', 'Calendar', 'Art', 'DateOfArt', 'DateOfHIVDiagnosis',
    'DateOfAIDSDiagnosis', 'DateOfArrival', 'DateOfBirth', 'AtRiskDate', 'U', 'Mig', 'KnownPrePost',
    'DTime', 'UniqueId', 'Consc', 'Consr', 'CobsTime', 'RobsTime', 'RLogObsTime2', 'Only'
  )
)
currentLevels <- levels(baseCD4VL2$Gender)
newLevels <- c('M', 'F')
levels(baseCD4VL2$Gender) <- newLevels[match(currentLevels, c('Male', 'Female'))]
baseCD4VL2[, Ord := data.table::rowid(UniqueId)]
baseCD4VL2[, UniqueId := UniqueId + max(baseAIDS$UniqueId)]
input <- list(
  Data = list(
    AIDS = baseAIDS,
    CD4VL = baseCD4VL2
  )
)

#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=icc
CCC=icpc
CXX=icpc
FC=ifort
AS=as

# Macros
CND_PLATFORM=Intel-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/computeFluxes.o \
	${OBJECTDIR}/generatePieceParDistri.o \
	${OBJECTDIR}/hadvppm.o \
	${OBJECTDIR}/secondOrderPolynomial.o \
	${OBJECTDIR}/unit-test-hadvppm.o \
	${OBJECTDIR}/updateConcentrations.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/pmcamx-hpc-sebastian

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/pmcamx-hpc-sebastian: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.f} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/pmcamx-hpc-sebastian ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/computeFluxes.o: computeFluxes.f03 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/computeFluxes.o computeFluxes.f03

${OBJECTDIR}/generatePieceParDistri.o: generatePieceParDistri.f03 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/generatePieceParDistri.o generatePieceParDistri.f03

${OBJECTDIR}/hadvppm.o: hadvppm.f03 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/hadvppm.o hadvppm.f03

${OBJECTDIR}/secondOrderPolynomial.o: secondOrderPolynomial.f03 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/secondOrderPolynomial.o secondOrderPolynomial.f03

${OBJECTDIR}/unit-test-hadvppm.o: unit-test-hadvppm.f03 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/unit-test-hadvppm.o unit-test-hadvppm.f03

${OBJECTDIR}/updateConcentrations.o: updateConcentrations.f03 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/updateConcentrations.o updateConcentrations.f03

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/pmcamx-hpc-sebastian
	${RM} *.mod

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc

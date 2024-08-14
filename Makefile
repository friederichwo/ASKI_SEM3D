#----------------------------------------------------------------------------
#   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
#   Copyright 2023 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
#
#   This file is part of ASKI version 1.2.
#
#   ASKI version 1.2 is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.
#
#   ASKI version 1.2 is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with ASKI version 1.2.  If not, see <http://www.gnu.org/licenses/>.
#----------------------------------------------------------------------------
#
################################################################
#  This is the Makefile for the ASKI main package (for GNU Make)
################################################################
#
#-----------------------------------------------------------------------
#  set the compiler
#
COMPILER = gfortran
CUDACOMPILER = /opt/nvidia/hpc_sdk/Linux_x86_64/22.11/compilers/bin/nvfortran
MPICOMPILER = mpif90
f2py_COMPILER = f2py
#
#-----------------------------------------------------------------------
#  General definitions
#
bindir = ./bin
obsdir = ./obj
#
FFLAGS = -O3 -J$(obsdir) -I/usr/include -Wunused-variable -Wuninitialized -fimplicit-none -ffixed-line-length-132 -fbounds-check -fbacktrace
#
#-----------------------------------------------------------------------
#  Direcories where to search for files to compile to .o by implicit rules below, and dependencies defined in rules.mk
#
vpath %.o $(obsdir)
vpath %.f90 ./f90
vpath %.cuf ./f90
#
#-----------------------------------------------------------------------
#  Implicit rule to compile .o files from .f90 files.
#  Because of vpath, targets and dependencies need not be
#  in the current directory.
#
%.o: %.f90
	$(COMPILER) -c $(FFLAGS) -I$(hdfinc) $< -o $(obsdir)/$@
%.o: %.cuf
	$(CUDACOMPILER) -c -Mpreprocess $< -o $(obsdir)/$@
#
#-----------------------------------------------------------------------
#  Object string for linking:
#  Adds object dir as prefix and removes directory part
#  of $^ (all dependencies)
#
obstring = $(addprefix $(obsdir)/,$(notdir $^))
#
#-----------------------------------------------------------------------
#  Library paths
#
# libraries for all applications: 
BLAS = /usr/lib/x86_64-linux-gnu/blas/libblas.so
BLAS_F2PY = -L/usr/libx86_64-linux-gnu -lblas
LAPACK = /usr/lib/x86_64-linux-gnu/liblapack.so
#
# libraries for parallel applications only:
BLACS = /usr/lib/x86_64-linux-gnu/libscalapack-openmpi.so
SCALAPACK = /usr/lib/x86_64-linux-gnu/libscalapack-openmpi.so
MPILIB = /usr/lib/x86_64-linux-gnu/openmpi/lib
#
#hdfinc = /usr/include
#hdflib = -L/usr/lib -lhdf5.dll
hdfinc = /usr/lib/x86_64-linux-gnu/hdf5/openmpi/include
hdflib = -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi/lib -lhdf5_fortran -lhdf5
hdflink = -Wl,-rpath,/usr/lib/x86_64-linux-gnu/hdf5/openmpi/lib
#
#  libraries for plotting
pgplot = -L/usr/lib -lpgplot $(x11)
x11 = -L/usr/lib/x86_64-linux-gnu -lX11

#
#-------------------------------------------------------------
#
.PHONY:
#
#----------------------------------------------------------------
#  Include dependencies:
#  dependencies is a Makefile because it is included. It containes all dependencies of
#  the .o files. If you change any such dependencies (e.g. by using an additional module
#  in some program/module), please update this file accordingly.
#
-include dependencies
#
#---------------------------------------------------------------
#
clean:
	-rm -f $(bindir)/*
	-rm -f $(obsdir)/*
	-rm -f ./mod/*
#
all: createCheckerboardKim computeKernelsDmspaceParallel computeMisfits computeModelNorm computeRegularSphericalGridSlices \
     collectColsums convertFmtomoToKim convertFmtomoToRsg decomposeDmspace evaluateSpmOnRsg initBasics kdispl2vtk kernel2vtk kgt2vtk \
     plotMeasuredWithSyntheticSeismograms solveCglsKernelSystem throwKimOnSphericalPseudoMesh \
     throwKernelOnSphericalPseudoMesh throwKimOnVgrid throwWavefieldOnSphericalPseudoMesh transformMeasuredDataAsciiSeis \
     transformSpecfemSynthetics updateModelPerturbations writeVtkFiles
#----------------------------------------------------------------
# Rules for all ASKI executables:
#
mpiSupport.o: mpiSupport.f90
	$(MPICOMPILER) -c -fintrinsic-modules-path $(MPILIB) -J$(obsdir) $< -o $(obsdir)/$@
hdfWrapper.o: hdfWrapper.f90
	$(MPICOMPILER) -c -I$(hdfinc) $(FFLAGS) $< -o $(obsdir)/$@
kernelLinearSystem.o: kernelLinearSystem.f90
	$(MPICOMPILER) -c -fintrinsic-modules-path $(MPILIB) -J$(obsdir) $< -o $(obsdir)/$@
linearModelRegularization.o: linearModelRegularization.f90
	$(MPICOMPILER) -c -fintrinsic-modules-path $(MPILIB) -J$(obsdir) $< -o $(obsdir)/$@
#
createCheckerboardKim: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o argumentParser.o \
        askiBackgroundModel.o axesRotation.o constants.o dateTime.o errorMessage.o fmtomoModel.o globalHdfInfo.o \
        globalMpiInfo.o hdfWrapper.o heapSort.o inputParameter.o inversionBasics.o inversionGrid.o \
        kernelInvertedModel.o locatePoint.o mathConstants.o propertySet.o readEventStationFile.o realloc.o \
        regularSphericalGrid.o seismicEvent.o seismicEventList.o seismicNetwork.o seismicStation.o \
        specfem3dForASKIFiles.o specfem3dInversionGrid.o string.o timeUtils.o
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#
computeHorizontalModelSlices: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o \
        argumentParser.o askiBackgroundModel.o axesRotation.o constants.o dateTime.o errorMessage.o fmtomoModel.o \
        globalHdfInfo.o globalMpiInfo.o hdfWrapper.o inputParameter.o inversionBasics.o inversionGrid.o \
        kernelInvertedModel.o mathConstants.o mpiSupport.o propertySet.o readEventStationFile.o realloc.o \
        seismicEvent.o seismicEventList.o seismicNetwork.o seismicStation.o smartUtils.o string.o timeUtils.o
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#
computeKernelsDmspaceParallel: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o \
        argumentParser.o axesRotation.o constants.o dataModelSpaceInfo.o dateTime.o errorMessage.o fileUnitHandler.o \
        globalHdfInfo.o globalMpiInfo.o heapSort.o hdfWrapper.o inputParameter.o inversionBasics.o inversionGrid.o \
        iterationStepBasics.o kernelWavefield.o locatePoint.o mathConstants.o mpiSupport.o propertySet.o rayInversionGrid.o \
        rayKernelWavefield.o readEventStationFile.o realloc.o seismicEvent.o seismicEventList.o seismicNetwork.o \
        seismicStation.o smartUtils.o specfem3dForASKIFiles.o specfem3dInversionGrid.o \
        specfem3dKernelWavefield.o spectralWaveformKernel.o string.o timeUtils.o
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#
computeMisfits: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o argumentParser.o \
	askiBackgroundModel.o axesRotation.o constants.o dataModelSpaceInfo.o dateTime.o errorMessage.o \
	fileUnitHandler.o fmtomoModel.o globalHdfInfo.o globalMpiInfo.o hdfWrapper.o inputParameter.o \
	inversionBasics.o inversionGrid.o kernelInvertedModel.o kernelLinearSystem.o kernelWavefield.o \
	mathConstants.o mpiSupport.o propertySet.o pythag.o rayKernelWavefield.o readEventStationFile.o realloc.o \
	seismicEvent.o seismicEventList.o seismicNetwork.o seismicStation.o specfem3dForASKIFiles.o \
	specfem3dKernelWavefield.o spectralWaveformKernel.o string.o timeUtils.o 
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink) $(BLAS) $(BLACS)
#
computeModelNorm: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o argumentParser.o \
        askiBackgroundModel.o axesRotation.o constants.o dataModelSpaceInfo.o dateTime.o errorMessage.o \
        fileUnitHandler.o fmtomoModel.o globalHdfInfo.o globalMpiInfo.o hdfWrapper.o heapSort.o inputParameter.o \
        inversionBasics.o inversionGrid.o iterationStepBasics.o kernelInvertedModel.o \
        linearModelRegularization.o locatePoint.o mathConstants.o mpiSupport.o propertySet.o \
        readEventStationFile.o realloc.o seismicEvent.o seismicEventList.o seismicNetwork.o seismicStation.o \
        smartUtils.o specfem3dForASKIFiles.o specfem3dInversionGrid.o string.o timeUtils.o
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink) $(BLAS) $(BLACS)
#
computeRegularSphericalGridSlices: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o \
	argumentParser.o axesRotation.o constants.o errorMessage.o globalHdfInfo.o hdfWrapper.o inputParameter.o \
	mathConstants.o realloc.o regularSphericalGrid.o string.o
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#
convertFmtomoToKim: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o argumentParser.o \
        askiBackgroundModel.o axesRotation.o constants.o dateTime.o errorMessage.o fileUnitHandler.o \
        fmtomoModel.o globalHdfInfo.o globalMpiInfo.o hdfWrapper.o heapSort.o inputParameter.o inversionBasics.o \
        inversionGrid.o iterationStepBasics.o kernelInvertedModel.o locatePoint.o mathConstants.o mpiSupport.o \
        propertySet.o readEventStationFile.o realloc.o seismicEvent.o seismicEventList.o seismicNetwork.o \
        seismicStation.o specfem3dForASKIFiles.o specfem3dInversionGrid.o string.o timeUtils.o
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#
convertFmtomoToRsg: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o argumentParser.o \
        axesRotation.o constants.o errorMessage.o fmtomoModel.o globalHdfInfo.o hdfWrapper.o mathConstants.o \
        realloc.o regularSphericalGrid.o string.o
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#
collectColsums: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o argumentParser.o \
        axesRotation.o constants.o dateTime.o errorMessage.o globalHdfInfo.o hdfWrapper.o inputParameter.o \
        inversionBasics.o mathConstants.o propertySet.o readEventStationFile.o realloc.o regularSphericalGrid.o \
        seismicEvent.o seismicEventList.o seismicNetwork.o seismicStation.o string.o timeUtils.o
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#
decomposeDmspace: %: %.o argumentParser.o constants.o dataModelSpaceInfo.o dateTime.o errorMessage.o \
	fileUnitHandler.o inputParameter.o inversionBasics.o iterationStepBasics.o mathConstants.o propertySet.o \
	readEventStationFile.o realloc.o seismicEvent.o seismicEventList.o seismicNetwork.o seismicStation.o \
	smartUtils.o string.o timeUtils.o 
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring)
#
evaluateSpmOnRsg: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o argumentParser.o \
        axesRotation.o constants.o errorMessage.o globalHdfInfo.o hdfWrapper.o mathConstants.o realloc.o \
        regularSphericalGrid.o string.o
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#
initBasics: %: %.o argumentParser.o constants.o dateTime.o errorMessage.o fileUnitHandler.o \
	inputParameter.o inversionBasics.o iterationStepBasics.o mathConstants.o propertySet.o \
	readEventStationFile.o realloc.o seismicEvent.o seismicEventList.o seismicNetwork.o seismicStation.o \
	string.o timeUtils.o 
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring)
#
kdispl2vtk: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o argumentParser.o axesRotation.o \
	constants.o dateTime.o errorMessage.o fileUnitHandler.o globalHdfInfo.o globalMpiInfo.o hdfWrapper.o \
	heapSort.o inputParameter.o inversionBasics.o inversionGrid.o iterationStepBasics.o kernelWavefield.o \
	locatePoint.o mathConstants.o mpiSupport.o propertySet.o readEventStationFile.o realloc.o seismicEvent.o \
	seismicEventList.o seismicNetwork.o seismicStation.o specfem3dForASKIFiles.o specfem3dInversionGrid.o \
	specfem3dKernelWavefield.o string.o timeUtils.o vtkFile.o 
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#
kernel2vtk: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o argumentParser.o axesRotation.o \
	constants.o dateTime.o errorMessage.o fileUnitHandler.o globalHdfInfo.o globalMpiInfo.o hdfWrapper.o \
	heapSort.o inputParameter.o inversionBasics.o inversionGrid.o iterationStepBasics.o kernelWavefield.o \
	locatePoint.o mathConstants.o mpiSupport.o propertySet.o rayKernelWavefield.o readEventStationFile.o \
	realloc.o seismicEvent.o seismicEventList.o seismicNetwork.o seismicStation.o specfem3dForASKIFiles.o \
	specfem3dInversionGrid.o specfem3dKernelWavefield.o spectralWaveformKernel.o string.o timeUtils.o \
	vtkFile.o 
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#
kgt2vtk: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o argumentParser.o axesRotation.o \
	constants.o dateTime.o errorMessage.o fileUnitHandler.o globalHdfInfo.o globalMpiInfo.o hdfWrapper.o \
	heapSort.o inputParameter.o inversionBasics.o inversionGrid.o iterationStepBasics.o kernelWavefield.o \
	locatePoint.o mathConstants.o mpiSupport.o propertySet.o readEventStationFile.o realloc.o seismicEvent.o \
	seismicEventList.o seismicNetwork.o seismicStation.o specfem3dForASKIFiles.o specfem3dInversionGrid.o \
	specfem3dKernelWavefield.o string.o timeUtils.o vtkFile.o 
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#
plotMeasuredWithSyntheticSeismograms: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o \
        argumentParser.o asciiDataIO.o asciiSynseisIO.o axesRotation.o basePlotGather.o changeAxesLimitsPgplot.o \
        constants.o dateTime.o errorMessage.o filterCoefficientsDecimate.o fourierSpectrum.o fourierTransform.o \
        hdfWrapper.o inputParameter.o interactiveBindingsPlotGather.o mathConstants.o pgPlotWindow.o \
        pythag.o readEventStationFile.o realloc.o recursiveFilterCoefficients.o seismicEvent.o seismicEventList.o \
        seismicNetwork.o seismicStation.o string.o timeSeries.o timeUtils.o
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(pgplot) $(hdflib) $(hdflink)
#
solveCglsKernelSystem: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o argumentParser.o \
	askiBackgroundModel.o axesRotation.o constants.o dataModelSpaceInfo.o dateTime.o errorMessage.o \
	fileUnitHandler.o fmtomoModel.o globalHdfInfo.o globalMpiInfo.o hdfWrapper.o heapSort.o inputParameter.o \
	inversionBasics.o inversionGrid.o iterationStepBasics.o kernelInvertedModel.o kernelLinearSystem.o \
	kernelWavefield.o linearModelRegularization.o locatePoint.o mathConstants.o mpiSupport.o propertySet.o \
	rayKernelWavefield.o readEventStationFile.o realloc.o seismicEvent.o seismicEventList.o seismicNetwork.o \
	pythag.o seismicStation.o smartUtils.o specfem3dForASKIFiles.o specfem3dInversionGrid.o \
	specfem3dKernelWavefield.o spectralWaveformKernel.o string.o timeUtils.o vtkFile.o 
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink) $(BLAS) $(BLACS)
#
throwKernelOnVgrid: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o argumentParser.o \
	axesRotation.o constants.o dateTime.o errorMessage.o fileUnitHandler.o globalHdfInfo.o globalMpiInfo.o \
	hdfWrapper.o heapSort.o inputParameter.o inversionBasics.o inversionGrid.o \
	kernelWavefield.o locatePoint.o mathConstants.o propertySet.o rayKernelWavefield.o \
	readEventStationFile.o realloc.o regularSphericalGrid.o seismicEvent.o seismicEventList.o \
	seismicNetwork.o seismicStation.o smartUtils.o specfem3dForASKIFiles.o specfem3dInversionGrid.o \
	specfem3dKernelWavefield.o spectralWaveformKernel.o string.o timeUtils.o 
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#
throwKernelOnSphericalPseudoMesh: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o \
	argumentParser.o axesRotation.o constants.o dateTime.o errorMessage.o fileUnitHandler.o globalHdfInfo.o \
	globalMpiInfo.o hdfWrapper.o heapSort.o inputParameter.o inversionBasics.o inversionGrid.o \
	iterationStepBasics.o kernelWavefield.o locatePoint.o mathConstants.o propertySet.o rayKernelWavefield.o \
	readEventStationFile.o realloc.o regularSphericalGrid.o seismicEvent.o seismicEventList.o \
	seismicNetwork.o seismicStation.o specfem3dForASKIFiles.o specfem3dInversionGrid.o \
	specfem3dKernelWavefield.o spectralWaveformKernel.o string.o timeUtils.o
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#
throwKimOnSphericalPseudoMesh: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o \
	argumentParser.o askiBackgroundModel.o axesRotation.o constants.o dateTime.o errorMessage.o \
	fileUnitHandler.o fmtomoModel.o globalHdfInfo.o globalMpiInfo.o hdfWrapper.o heapSort.o inputParameter.o \
	inversionBasics.o inversionGrid.o iterationStepBasics.o kernelInvertedModel.o locatePoint.o \
	mathConstants.o propertySet.o readEventStationFile.o realloc.o regularSphericalGrid.o seismicEvent.o \
	seismicEventList.o seismicNetwork.o seismicStation.o smartUtils.o specfem3dForASKIFiles.o \
	specfem3dInversionGrid.o string.o timeUtils.o
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#
throwKimOnVgrid: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o argumentParser.o \
	askiBackgroundModel.o axesRotation.o constants.o dateTime.o errorMessage.o fileUnitHandler.o \
	fmtomoModel.o globalHdfInfo.o globalMpiInfo.o hdfWrapper.o heapSort.o inputParameter.o inversionBasics.o \
	inversionGrid.o kernelInvertedModel.o locatePoint.o mathConstants.o propertySet.o \
	readEventStationFile.o realloc.o regularSphericalGrid.o seismicEvent.o seismicEventList.o \
	seismicNetwork.o seismicStation.o smartUtils.o specfem3dForASKIFiles.o specfem3dInversionGrid.o string.o \
	timeUtils.o 
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#
throwWavefieldOnSphericalPseudoMesh: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o \
        argumentParser.o axesRotation.o constants.o dateTime.o errorMessage.o fileUnitHandler.o globalHdfInfo.o \
        globalMpiInfo.o hdfWrapper.o heapSort.o inputParameter.o inversionBasics.o inversionGrid.o \
        iterationStepBasics.o kernelWavefield.o locatePoint.o mathConstants.o mpiSupport.o propertySet.o rayKernelWavefield.o \
        readEventStationFile.o realloc.o regularSphericalGrid.o seismicEvent.o seismicEventList.o \
        seismicNetwork.o seismicStation.o specfem3dForASKIFiles.o specfem3dInversionGrid.o \
        specfem3dKernelWavefield.o string.o timeUtils.o
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#
transformMeasuredDataAsciiSeis: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o \
        argumentParser.o asciiSynseisIO.o axesRotation.o constants.o dateTime.o discreteFourierTransform.o errorMessage.o \
        hdfWrapper.o inputParameter.o inversionBasics.o locatePoint.o mathConstants.o propertySet.o readEventStationFile.o \
        realloc.o seismicEvent.o seismicEventList.o seismicNetwork.o seismicStation.o smartUtils.o \
        string.o timeUtils.o travelTimes.o
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#
transformSpecfemSynthetics: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o \
        argumentParser.o asciiDataIO.o asciiSynseisIO.o axesRotation.o constants.o dateTime.o discreteFourierTransform.o \
        errorMessage.o hdfWrapper.o inputParameter.o inversionBasics.o locatePoint.o mathConstants.o propertySet.o \
        readEventStationFile.o realloc.o seismicEvent.o seismicEventList.o seismicNetwork.o seismicStation.o \
        smartUtils.o string.o timeUtils.o travelTimes.o
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#
updateModelPerturbations: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o argumentParser.o \
        askiBackgroundModel.o axesRotation.o constants.o dateTime.o errorMessage.o fmtomoModel.o globalHdfInfo.o \
        globalMpiInfo.o hdfWrapper.o heapSort.o inputParameter.o inversionBasics.o inversionGrid.o \
        kernelInvertedModel.o locatePoint.o mathConstants.o propertySet.o readEventStationFile.o realloc.o \
        seismicEvent.o seismicEventList.o seismicNetwork.o seismicStation.o specfem3dForASKIFiles.o \
        specfem3dInversionGrid.o string.o timeUtils.o
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#
writeVtkFiles: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o argumentParser.o \
	axesRotation.o constants.o dateTime.o errorMessage.o fileUnitHandler.o globalHdfInfo.o globalMpiInfo.o \
	hdfWrapper.o heapSort.o inputParameter.o inversionBasics.o inversionGrid.o iterationStepBasics.o \
	locatePoint.o mathConstants.o mpiSupport.o propertySet.o readEventStationFile.o realloc.o seismicEvent.o \
	seismicEventList.o seismicNetwork.o seismicStation.o specfem3dForASKIFiles.o specfem3dInversionGrid.o \
	string.o timeUtils.o vtkFile.o 
	$(MPICOMPILER) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#
testcudafor: %: %.o
	$(CUDACOMPILER) -o $(bindir)/$@ $(obstring) -L/usr/local/cuda-12/lib64 -lcudart
#
fpy: fortranForPython.f90
	$(f2py_COMPILER) -c -m fpy f90/fortranForPython.f90



    

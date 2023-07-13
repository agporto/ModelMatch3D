##BASE PYTHON
import os
import logging
import json
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import vtk.util.numpy_support as vtk_np
import numpy as np

#
# BiteID
#

class BiteID(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "BiteID"
    self.parent.categories = ["Examples"]
    self.parent.dependencies = []
    self.parent.contributors = ["Arthur Porto"]
    self.parent.helpText = """
     This module automatically registers two dental scans, establishing point correspondences between them and computing metrics of overall similarity for individual identification purposes.
See more information in <a href="https://github.com/organization/projectname#IdentiScan">module documentation</a>.
      """
    self.parent.acknowledgementText = """
This file was originally developed by Arthur Porto (add others) and it was originally supported by a XXXX grant.
      """
#
# BiteIDWidget
#

class BiteIDWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

  def __init__(self, parent=None):
    ScriptedLoadableModuleWidget.__init__(self, parent)

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)
    # Ensure that correct version of open3d Python package is installed
    needRestart = False
    needInstall = False
    Open3dVersion = "0.17" if slicer.app.os != "macosx" else "0.16.1"
    try:
      import open3d as o3d
      from packaging import version
      if version.parse(o3d.__version__) != version.parse(Open3dVersion):
        if not slicer.util.confirmOkCancelDisplay(f"BiteID requires installation of open3d (version {Open3dVersion}).\nClick OK to upgrade open3d and restart the application."):
          self.ui.showBrowserOnEnter = False
          return
        needRestart = True
        needInstall = True
    except ModuleNotFoundError:
      needInstall = True

    if needInstall:
      progressDialog = slicer.util.createProgressDialog(labelText='Upgrading open3d. This may take a minute...', maximum=0)
      slicer.app.processEvents()
      # the argument into multiple command-line arguments when there are spaces in the path)
      slicer.util.pip_install(f'open3d=={Open3dVersion}')
      import open3d as o3d
      progressDialog.close()
    if needRestart:
      slicer.util.restart()

    # Load widget from .ui file (created by Qt Designer).
    uiWidget = slicer.util.loadUI(self.resourcePath('UI/BiteID.ui'))
    self.layout.addWidget(uiWidget)
    self.ui = slicer.util.childWidgetVariables(uiWidget)

    # Setup checkbox stylesheets
    moduleDir = os.path.dirname(slicer.util.modulePath(self.__module__))
    self.onIconPath = moduleDir +'/Resources/Icons/switch_on_v2.png'
    self.offIconPath = moduleDir +'/Resources/Icons/switch_off_v2.png'

    # Single Alignment connections
    self.ui.sourceModelSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.ui.sourceModelSelector.setMRMLScene( slicer.mrmlScene)
    self.ui.targetModelSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.ui.targetModelSelector.setMRMLScene( slicer.mrmlScene)
    self.ui.runBiteIDButton.connect('clicked(bool)', self.onRunBiteIDButton) #Connection for run all BiteID steps button
    self.ui.switchSettingsButton.connect('clicked(bool)', self.onSwitchSettingsButton)
    self.ui.showTargetModelCheckBox.connect('toggled(bool)', self.onshowTargetModelCheckBox)
    self.setCheckboxStyle(self.ui.showTargetModelCheckBox)
    self.ui.showSourceModelCheckBox.connect('toggled(bool)', self.onshowSourceModelCheckBox)
    self.setCheckboxStyle(self.ui.showSourceModelCheckBox)
    self.ui.signedDistanceCheckBox.connect('toggled(bool)', self.onshowSignedDistanceMapCheckBox)
    self.setCheckboxStyle(self.ui.signedDistanceCheckBox)


    # Advanced Settings connections
    self.ui.pointDensityAdvancedSlider.connect('valueChanged(double)', self.onChangeAdvanced)
    self.ui.normalSearchRadiusSlider.connect('valueChanged(double)', self.onChangeAdvanced)
    self.ui.FPFHSearchRadiusSlider.connect('valueChanged(double)', self.onChangeAdvanced)
    self.ui.maximumCPDThreshold.connect('valueChanged(double)', self.onChangeAdvanced)
    self.ui.maxRANSAC.connect('valueChanged(double)', self.onChangeAdvanced)
    self.ui.RANSACConfidence.connect('valueChanged(double)', self.onChangeAdvanced)
    self.ui.ICPDistanceThresholdSlider.connect('valueChanged(double)', self.onChangeAdvanced)


    # Batch Processing connections
    self.ui.methodBatchWidget.connect('currentIndexChanged(int)', self.onChangeMultiTemplate)
    self.ui.sourceModelMultiSelector.connect('validInputChanged(bool)', self.onSelectMultiProcess)
    self.ui.targetModelMultiSelector.connect('validInputChanged(bool)', self.onSelectMultiProcess)
    self.ui.landmarkOutputSelector.connect('validInputChanged(bool)', self.onSelectMultiProcess)
    self.ui.applyLandmarkMultiButton.connect('clicked(bool)', self.onApplyLandmarkMulti)


    # initialize the parameter dictionary from single run parameters
    self.parameterDictionary = {
      "pointDensity": self.ui.pointDensityAdvancedSlider.value,
      "normalSearchRadius" : self.ui.normalSearchRadiusSlider.value,
      "FPFHSearchRadius" : self.ui.FPFHSearchRadiusSlider.value,
      "distanceThreshold" : self.ui.maximumCPDThreshold.value,
      "maxRANSAC" : int(self.ui.maxRANSAC.value),
      "RANSACConfidence" : self.ui.RANSACConfidence.value,
      "ICPDistanceThreshold"  : self.ui.ICPDistanceThresholdSlider.value,
      }
    #Set up table
    #self.tableNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTableNode')
    self.tableNode = slicer.vtkMRMLTableNode()
    # Add columns to the table
    self.tableNode.AddColumn().SetName('Query')
    self.tableNode.AddColumn().SetName('Reference')
    self.tableNode.AddColumn().SetName('Fitness')
    self.tableNode.AddColumn().SetName('Partial Fitness')
    self.tableNode.AddColumn().SetName('RMSE')
    self.tableNode.AddColumn().SetName('Partial RMSE')
    self.tableNode.SetUseColumnNameAsColumnHeader(True)

    self.ui.potentialMatches.setMRMLTableNode(self.tableNode)
    self.ui.potentialMatches.resizeColumnsToContents()
    self.ui.potentialMatches.horizontalHeader().setSectionResizeMode(qt.QHeaderView.Stretch)

    self.lineNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsLineNode', "L")
    
    self.ui.rulerWidget.placeButton().show()
    self.ui.rulerWidget.deleteButton().show()
    self.ui.rulerWidget.setDeleteAllControlPointsOptionVisible(True)
    self.ui.rulerWidget.setCurrentNode(self.lineNode)
    self.ui.rulerWidget.setMRMLScene(slicer.mrmlScene)
    self.lineNode.GetDisplayNode().SetVisibility(True)
    #self.ui.rulerWidget.show()



  def cleanup(self):
    pass

  def setCheckboxStyle(self, checkbox):
    checkbox.setStyleSheet("QCheckBox::indicator:unchecked{image: url(" + self.offIconPath + "); } \n"
      + "QCheckBox::indicator:checked{image: url(" + self.onIconPath + "); } \n"
      + "QCheckBox::indicator{width: 50px;height: 25px;}\n")

  def onSelect(self):
    self.ui.runBiteIDButton.enabled = bool ( self.ui.sourceModelSelector.currentNode() and self.ui.targetModelSelector.currentNode())
    self.ui.switchSettingsButton.enabled = bool ( self.ui.sourceModelSelector.currentNode() and self.ui.targetModelSelector.currentNode())

  def onSelectMultiProcess(self):
    self.ui.applyLandmarkMultiButton.enabled = bool ( self.ui.sourceModelMultiSelector.currentPath and self.ui.targetModelMultiSelector.currentPath and self.ui.landmarkOutputSelector.currentPath )
  

  def onRunBiteIDButton(self):

    logic = BiteIDLogic()
    logic.hideScene()
    self.sourceModelNode_orig = self.ui.sourceModelSelector.currentNode()
    self.sourceModelNode_orig.GetDisplayNode().SetVisibility(False)
    #Clone the original source mesh stored in the node sourceModelNode_orig
    shNode = slicer.vtkMRMLSubjectHierarchyNode.GetSubjectHierarchyNode(slicer.mrmlScene)
    itemIDToClone = shNode.GetItemByDataNode(self.sourceModelNode_orig)
    clonedItemID = slicer.modules.subjecthierarchy.logic().CloneSubjectHierarchyItem(shNode, itemIDToClone)
    self.sourceModelNode = shNode.GetItemDataNode(clonedItemID)
    self.sourceModelNode.GetDisplayNode().SetVisibility(True)
    self.sourceModelNode.SetName("Query model(registered)")#+run_counter) #Create a cloned source model node
    # slicer.mrmlScene.RemoveNode(self.sourceModelNode_copy)
    #
    self.targetModelNode = self.ui.targetModelSelector.currentNode()
    self.targetModelNode.GetDisplayNode().SetVisibility(True)
    self.updateLayout()
    #
    self.sourcePoints, self.targetPoints, self.sourceFeatures, \
      self.targetFeatures, self.voxelSize = logic.runSubsample(self.sourceModelNode, self.targetModelNode, self.parameterDictionary)
    
    # Convert to VTK points for visualization
    self.targetVTK = logic.convertPointsToVTK(self.targetPoints.points)

    blue=[0,0,1]
    self.targetCloudNode = logic.displayPointCloud(self.targetVTK, self.voxelSize / 5, 'Target Pointcloud', blue)
    self.targetCloudNode.GetDisplayNode().SetVisibility(False)

    #
    #RANSAC & ICP transformation of source pointcloud
    self.transformMatrix, self.fitness, self.partialFitness, self.inlier_rmse, self.partial_inlier_rmse = logic.estimateTransform(self.sourcePoints, self.targetPoints, self.sourceFeatures, self.targetFeatures, self.voxelSize, self.parameterDictionary)
    self.ui.totalFitness.value = self.fitness*100
    self.ui.partialFitness.value = self.partialFitness*100

    # Custom color for progress bar
    color_pallete = ['#3cae3c', '#e6e600', '#f08080']#['#d2f8d2', '#e6e600', '#f08080']
    color_total = color_pallete[0] if self.fitness > 0.3 else color_pallete[1] if self.fitness > 0.15 else color_pallete[2]
    color_partial = color_pallete[0] if self.partialFitness > 0.8 else color_pallete[1] if self.partialFitness > 0.4 else color_pallete[2]
    self.ui.totalFitness.setStyleSheet(f"""
      QProgressBar{{
          border: 1px solid grey;
          border-radius: 5px;
          text-align: center
      }}

      QProgressBar::chunk {{
          background-color: {color_total};
          width: 10px;
          margin: 0px;
      }}
      """)
    self.ui.partialFitness.setStyleSheet(f"""
      QProgressBar{{
          border: 1px solid grey;
          border-radius: 5px;
          text-align: center
      }}

      QProgressBar::chunk {{
          background-color: {color_partial};
          width: 10px;
          margin: 0px;
      }}
      """)
    #
    self.ICPTransformNode = logic.convertMatrixToTransformNode(self.transformMatrix, 'Rigid Transformation Matrix')#+ run_counter)
    self.sourcePoints.transform(self.transformMatrix)

    #Setup source pointcloud VTK object
    self.sourceVTK = logic.convertPointsToVTK(self.sourcePoints.points)
    #
    red=[1,0,0]
    self.sourceCloudNode = logic.displayPointCloud(self.sourceVTK, self.voxelSize / 5, ('Query Pointcloud (registered)'), red)#+ run_counter), red)
    self.sourceCloudNode.GetDisplayNode().SetVisibility(False)
    #
    #Transform the source model
    self.sourceModelNode.SetAndObserveTransformNodeID(self.ICPTransformNode.GetID())
    slicer.vtkSlicerTransformLogic().hardenTransform(self.sourceModelNode)
    #
    self.sourceModelNode.GetDisplayNode().SetColor(red)

    #Create a folder and add all objects into the folder
    folderNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    sceneItemID = folderNode.GetSceneItemID()
    folderName = 'BiteID_output'
    existingFolder = folderNode.GetItemByName(folderName)
    if existingFolder:
      folderNode.RemoveItem(existingFolder)

    newFolder = folderNode.CreateFolderItem(sceneItemID , folderName)#+ run_counter)
    sourceModelItem = folderNode.GetItemByDataNode(self.sourceModelNode)
    sourceCloudItem = folderNode.GetItemByDataNode(self.sourceCloudNode)
    targetCloudItem = folderNode.GetItemByDataNode(self.targetCloudNode)
    ICPTransformItem = folderNode.GetItemByDataNode(self.ICPTransformNode)

    
    folderNode.SetItemParent(sourceModelItem, newFolder)
    folderNode.SetItemParent(sourceCloudItem, newFolder)
    folderNode.SetItemParent(targetCloudItem, newFolder)
    folderNode.SetItemParent(ICPTransformItem, newFolder)

    #
    #Enable visualization
    #
    self.ui.showTargetModelCheckBox.enabled = True
    self.ui.showTargetModelCheckBox.checked = 1
    #
    self.ui.showSourceModelCheckBox.enabled = True
    self.ui.showSourceModelCheckBox.checked = 1

    self.ui.signedDistanceCheckBox.enabled = True
    self.ui.signedDistanceCheckBox.checked = 0
    self.lineNode.GetDisplayNode().SetVisibility(True)




  def onSwitchSettingsButton(self):
    self.ui.tabsWidget.setCurrentWidget(self.ui.advancedSettingsTab)
    self.ui.runBiteIDButton.enabled = True

  def onshowTargetModelCheckBox(self):
    try:
      if self.ui.showTargetModelCheckBox.isChecked():
        self.targetModelNode.GetDisplayNode().SetVisibility(True)
      else:
        self.targetModelNode.GetDisplayNode().SetVisibility(False)
    except:
      self.ui.showTargetModelCheckBox.enabled = False

  def onshowSourceModelCheckBox(self):
    try:
      if self.ui.showSourceModelCheckBox.isChecked():
        self.sourceModelNode.GetDisplayNode().SetVisibility(True)
      else:
        self.sourceModelNode.GetDisplayNode().SetVisibility(False)
    except:
      self.ui.showSourceModelCheckBox.enabled = False
  
  def onshowSignedDistanceMapCheckBox(self):
    logic = BiteIDLogic()
    try:
      if self.ui.signedDistanceCheckBox.isChecked():
        self.targetModelNode.GetDisplayNode().SetVisibility(False)
        self.sourceModelNode = logic.signedDistancePainting(self.sourceModelNode, self.targetModelNode, self.voxelSize)
      else:
        self.sourceModelNode.GetDisplayNode().SetVisibility(True)
        self.sourceModelNode.GetDisplayNode().SetScalarVisibility(False)
        red = [1,0,0]
        self.sourceModelNode.GetDisplayNode().SetColor(red)
        self.targetModelNode.GetDisplayNode().SetVisibility(True)
    except:
      pass

  def onApplyLandmarkMulti(self):
    table = self.tableNode.GetTable()
    table.SetNumberOfRows(0)
    self.ui.potentialMatches.update()
    slicer.app.processEvents()

    logic = BiteIDLogic()
    sourceModelPath = self.ui.sourceModelMultiSelector.currentPath
    targetModelDirectory = self.ui.targetModelMultiSelector.currentPath
    outputDirectory = self.ui.landmarkOutputSelector.currentPath

    extension= ".csv"
    TargetModelList = []
    extras = {"Source" : sourceModelPath, "Target" : targetModelDirectory, "Output" : outputDirectory}
    extras.update(self.parameterDictionary)
    parameterFile = os.path.join(outputDirectory, 'queryParameters.txt')
    json.dump(extras, open(parameterFile,'w'), indent = 2)

    for targetFileName in os.listdir(targetModelDirectory):
        if targetFileName.endswith((".ply", ".obj", ".vtk", ".stl")):
            targetFilePath = os.path.join(targetModelDirectory, targetFileName)
            TargetModelList.append(targetFilePath)

    allComparisonsDir = os.path.join(outputDirectory, 'allComparisons')
    os.makedirs(allComparisonsDir, exist_ok=True)

    if os.path.isdir(sourceModelPath):
        sourceFiles = [os.path.join(sourceModelPath, f) for f in os.listdir(sourceModelPath) if f.endswith((".ply", ".obj", ".vtk", ".stl"))]
    elif os.path.isfile(sourceModelPath):
        sourceFiles = [sourceModelPath]
    else:
        print('::::Could not find the file or directory in question')
        return
    totalComparisons = len(sourceFiles) * len(TargetModelList)
    count = 0
    for sourceFilePath in sourceFiles:
        basename = os.path.splitext(os.path.basename(sourceFilePath))[0]
        allComparisons = []
        
        for targetFilePath in TargetModelList:
            count += 1
            referenceName = os.path.splitext(os.path.basename(targetFilePath))[0]
            resultList = logic.pairwiseAlignment(sourceFilePath, targetFilePath, self.parameterDictionary)
            resultList = [ round(elem, 2) for elem in resultList ]

            # Append the results to allComparisons list along with target model path
            allComparisons.append((referenceName, resultList))


            # print highly similar models to the Slicer GUI
            if resultList[0] > 0.17 and resultList[1] > 0.55:  # adjust the threshold as needed
              newRowId = table.InsertNextBlankRow() 
              data = [basename, referenceName] + resultList
              print(data)
              for col in range(len(data)):
                table.SetValue(newRowId, col, data[col])
              table.Modified()
              self.ui.potentialMatches.update()

            self.ui.progressBar.setValue(100 * count / totalComparisons)
            slicer.app.processEvents()



        # save all comparisons to a single report file
        reportFilePath = os.path.join(allComparisonsDir, basename + extension)
        with open(reportFilePath, 'w') as reportFile:
            reportFile.write('Reference Model,Fitness,Partial Fitness,RMSE,Partial RMSE\n')
            for modelPath, results in sorted(allComparisons, key=lambda x: x[1][0], reverse=True):
                reportFile.write(f'{modelPath},{",".join(map(str, results))}\n')


  def updateLayout(self):
    layoutManager = slicer.app.layoutManager()
    layoutManager.setLayout(9)  #set layout to 3D only
    layoutManager.threeDWidget(0).threeDView().resetFocalPoint()
    layoutManager.threeDWidget(0).threeDView().resetCamera()

  def onChangeAdvanced(self):
    self.updateParameterDictionary()

  def onChangeMultiTemplate(self):
    if self.ui.methodBatchWidget.currentText == 'Multi-Query':
      self.ui.sourceModelMultiSelector.currentPath = ''
      self.ui.sourceModelMultiSelector.filters  = ctk.ctkPathLineEdit.Dirs
      self.ui.sourceModelMultiSelector.toolTip = "Select a directory containing the query models"

    else:
      self.ui.sourceModelMultiSelector.filters  = ctk.ctkPathLineEdit().Files
      self.ui.sourceModelMultiSelector.nameFilters  = ["*.ply *.obj *.vtk *.stl"]
      self.ui.sourceModelMultiSelector.toolTip = "Select the query model"


  def updateParameterDictionary(self):
    if hasattr(self, 'parameterDictionary'):
      self.parameterDictionary["pointDensity"] = self.ui.pointDensityAdvancedSlider.value
      self.parameterDictionary["normalSearchRadius"] = int(self.ui.normalSearchRadiusSlider.value)
      self.parameterDictionary["FPFHSearchRadius"] = int(self.ui.FPFHSearchRadiusSlider.value)
      self.parameterDictionary["distanceThreshold"] = self.ui.maximumCPDThreshold.value
      self.parameterDictionary["maxRANSAC"] = int(self.ui.maxRANSAC.value)
      self.parameterDictionary["RANSACConfidence"] = self.ui.RANSACConfidence.value
      self.parameterDictionary["ICPDistanceThreshold"] = self.ui.ICPDistanceThresholdSlider.value

#
# BiteIDLogic
#

class BiteIDLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
    computation done by your module.  The interface
    should be such that other python code can import
    this class and make use of the functionality without
    requiring an instance of the Widget.
    Uses ScriptedLoadableModuleLogic base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """
  def pairwiseAlignment (self, sourceFilePath, targetFilePath, parameters):
    targetModelNode = slicer.util.loadModel(targetFilePath)
    targetModelNode.GetDisplayNode().SetVisibility(False)
    sourceModelNode = slicer.util.loadModel(sourceFilePath)
    sourceModelNode.GetDisplayNode().SetVisibility(False)
    sourcePoints, targetPoints, sourceFeatures, targetFeatures, voxelSize = self.runSubsample(sourceModelNode,
        targetModelNode, parameters)
    _, icpFitness, ransacFitness, icpRmse, ransacRmse = self.estimateTransform(sourcePoints, targetPoints, sourceFeatures, targetFeatures, voxelSize, parameters)
    slicer.mrmlScene.RemoveNode(sourceModelNode)
    slicer.mrmlScene.RemoveNode(targetModelNode)
    return [icpFitness, ransacFitness, icpRmse, ransacRmse]

  def hideScene(self):
    allNodes = slicer.mrmlScene.GetNodes()
    for i in range(allNodes.GetNumberOfItems()):
        node = allNodes.GetItemAsObject(i)
        if node.IsA('vtkMRMLDisplayableNode'):
            node.SetDisplayVisibility(False)

  def convertMatrixToVTK(self, matrix):
    matrix_vtk = vtk.vtkMatrix4x4()
    for i in range(4):
      for j in range(4):
        matrix_vtk.SetElement(i,j,matrix[i][j])
    return matrix_vtk

  def convertMatrixToTransformNode(self, matrix, transformName):
    matrix_vtk = vtk.vtkMatrix4x4()
    for i in range(4):
      for j in range(4):
        matrix_vtk.SetElement(i,j,matrix[i][j])

    transform = vtk.vtkTransform()
    transform.SetMatrix(matrix_vtk)
    transformNode =  slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTransformNode', transformName)
    transformNode.SetAndObserveTransformToParent( transform )

    return transformNode

  def applyTransform(self, matrix, polydata):
    transform = vtk.vtkTransform()
    transform.SetMatrix(matrix)

    transformFilter = vtk.vtkTransformPolyDataFilter()
    transformFilter.SetTransform(transform)
    transformFilter.SetInputData(polydata)
    transformFilter.Update()
    return transformFilter.GetOutput()

  def convertPointsToVTK(self, points):
    array_vtk = vtk_np.numpy_to_vtk(points, deep=True, array_type=vtk.VTK_FLOAT)
    points_vtk = vtk.vtkPoints()
    points_vtk.SetData(array_vtk)
    polydata_vtk = vtk.vtkPolyData()
    polydata_vtk.SetPoints(points_vtk)
    return polydata_vtk


  def displayPointCloud(self, polydata, pointRadius, nodeName, nodeColor):
    #set up glyph for visualizing point cloud
    sphereSource = vtk.vtkSphereSource()
    sphereSource.SetRadius(pointRadius)
    glyph = vtk.vtkGlyph3D()
    glyph.SetSourceConnection(sphereSource.GetOutputPort())
    glyph.SetInputData(polydata)
    glyph.ScalingOff()
    glyph.Update()

    #display
    modelNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode', nodeName)
    modelNode.CreateDefaultDisplayNodes()

    modelNode.SetAndObservePolyData(glyph.GetOutput())
    modelNode.GetDisplayNode().SetColor(nodeColor)
    return modelNode
  
  def signedDistancePainting(self, sourceModel, targetModel, voxelSize):
    #get distance
    distance_filter = vtk.vtkDistancePolyDataFilter()
    distance_filter.SetInputData(0, sourceModel.GetPolyData())
    distance_filter.SetInputData(1, targetModel.GetPolyData())
    distance_filter.Update()


    # Get distance values from the filter output
    distance_data = distance_filter.GetOutput().GetPointData().GetScalars()

    # Convert distance data to a numpy array
    distance_array = vtk_np.vtk_to_numpy(distance_data)

    # Calculate mean and standard deviation
    mean_distance = np.mean(distance_array)
    std_distance = np.std(distance_array)

    # Print the mean and standard deviation
    print("Mean distance:", mean_distance)
    print("Standard deviation:", std_distance)



    #display
    sourceModel.SetAndObservePolyData(distance_filter.GetOutput())
    sourceModel.GetDisplayNode().SetActiveScalarName('Distance')
    sourceModel.GetDisplayNode().SetAndObserveColorNodeID(slicer.util.getNode('fMRIPA').GetID())
    sourceModel.GetDisplayNode().SetScalarVisibility(True)

    # Set scalar range
    sourceModel.GetDisplayNode().SetScalarRangeFlag(0)
    sourceModel.GetDisplayNode().SetScalarRange(-voxelSize, voxelSize)

    return sourceModel
  
  def signedDistancePainting2(self, sourceModel, targetModel, voxelSize):
      from open3d import geometry, utility

      # Convert source and target models to Open3D point clouds
      sourcePoints = slicer.util.arrayFromModelPoints(sourceModel)
      source = geometry.PointCloud()
      source.points = utility.Vector3dVector(sourcePoints)
      targetPoints = slicer.util.arrayFromModelPoints(targetModel)
      target = geometry.PointCloud()
      target.points = utility.Vector3dVector(targetPoints)

      # Compute point cloud distance using Open3D
      distance = source.compute_point_cloud_distance(target)
      distance = np.asarray(distance)

      mean = np.mean(distance)
      distance = distance - mean
      std = np.std(distance)
      print("Mean: " + str(mean))
      print("Std: " + str(std))


      # Create a vtkDoubleArray and set the scalar values
      distanceArray = vtk.vtkDoubleArray()
      distanceArray.SetName("Distance")
      distanceArray.SetNumberOfComponents(1)
      distanceArray.SetArray(distance, distance.size, 1)

      # Set the scalar array to the point data of the source model
      sourceModel.GetPolyData().GetPointData().SetScalars(distanceArray)

      # Set up the display properties for the source model node
      sourceModel.GetDisplayNode().SetActiveScalarName("Distance")
      sourceModel.GetDisplayNode().SetAndObserveColorNodeID(slicer.util.getNode("fMRIPA").GetID())
      sourceModel.GetDisplayNode().SetScalarVisibility(True)

      # Set scalar range
      sourceModel.GetDisplayNode().SetScalarRangeFlag(0)
      sourceModel.GetDisplayNode().SetScalarRange(-voxelSize, voxelSize)

      return sourceModel


  
  def displayMesh(self, polydata, nodeName, nodeColor):
    modelNode=slicer.mrmlScene.GetFirstNodeByName(nodeName)
    if modelNode is None:  # if there is no node with this name, create with display node
      modelNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode', nodeName)
      modelNode.CreateDefaultDisplayNodes()

    modelNode.SetAndObservePolyData(polydata)
    modelNode.GetDisplayNode().SetColor(nodeColor)
    return modelNode

  def estimateTransform(self, sourcePoints, targetPoints, sourceFeatures, targetFeatures, voxelSize, parameters):
    ransac = self.execute_global_registration(sourcePoints, targetPoints, sourceFeatures, targetFeatures, voxelSize,
      parameters["distanceThreshold"], parameters["maxRANSAC"], parameters["RANSACConfidence"])
    # Refine the initial registration using an Iterative Closest Point (ICP) registration
    #import time
    icp = self.refine_registration(sourcePoints, targetPoints, voxelSize, ransac, parameters["ICPDistanceThreshold"])
    return icp.transformation, icp.fitness, ransac.fitness, icp.inlier_rmse, ransac.inlier_rmse

  def runSubsample(self, sourceModel, targetModel, parameters):
    from open3d import io
    from open3d import geometry
    from open3d import utility
    print(":: Loading point clouds and downsampling")
    sourcePoints = slicer.util.arrayFromModelPoints(sourceModel)
    source = geometry.PointCloud()
    source.points = utility.Vector3dVector(sourcePoints)
    targetPoints = slicer.util.arrayFromModelPoints(targetModel)
    target = geometry.PointCloud()
    target.points = utility.Vector3dVector(targetPoints)
    targetSize = np.linalg.norm(np.asarray(target.get_max_bound()) - np.asarray(target.get_min_bound()))
    voxel_size = targetSize/(130*parameters["pointDensity"])
    source_down, source_fpfh = self.preprocess_point_cloud(source, voxel_size, parameters["normalSearchRadius"], parameters["FPFHSearchRadius"])
    target_down, target_fpfh = self.preprocess_point_cloud(target, voxel_size, parameters["normalSearchRadius"], parameters["FPFHSearchRadius"])
    return source_down, target_down, source_fpfh, target_fpfh, voxel_size


  def preprocess_point_cloud(self, pcd, voxel_size, radius_normal_factor, radius_feature_factor):
    from open3d import geometry
    from open3d import pipelines
    registration = pipelines.registration
    print(":: Downsample with a voxel size %.3f." % voxel_size)
    pcd_down = pcd.voxel_down_sample(voxel_size)
    radius_normal = voxel_size * radius_normal_factor
    print(":: Estimate normal with search radius %.3f." % radius_normal)
    pcd_down.estimate_normals(
        geometry.KDTreeSearchParamHybrid(radius=radius_normal, max_nn=30))
    radius_feature = voxel_size * radius_feature_factor
    print(":: Compute FPFH feature with search radius %.3f." % radius_feature)
    pcd_fpfh = registration.compute_fpfh_feature(
        pcd_down,
        geometry.KDTreeSearchParamHybrid(radius=radius_feature, max_nn=100))
    return pcd_down, pcd_fpfh


  def execute_global_registration(self, source_down, target_down, source_fpfh,
                                target_fpfh, voxel_size, distance_threshold_factor, maxIter, confidence):
    from open3d import pipelines
    from open3d import utility
    utility.random.seed(99)
    registration = pipelines.registration
    distance_threshold = voxel_size * distance_threshold_factor
    print(":: RANSAC registration on downsampled point clouds.")
    print("   Since the downsampling voxel size is %.3f," % voxel_size)
    print("   we use a liberal distance threshold %.3f." % distance_threshold)
    fitness = 0
    count = 0
    maxAttempts = 3
    while fitness < 0.99 and count < maxAttempts:
        result = registration.registration_ransac_based_on_feature_matching(
        source_down, target_down, source_fpfh, target_fpfh, True,
        distance_threshold,
       registration.TransformationEstimationPointToPoint(False),
        3, [
            registration.CorrespondenceCheckerBasedOnEdgeLength(
                0.9 ),
            registration.CorrespondenceCheckerBasedOnDistance(
                distance_threshold)
        ], registration.RANSACConvergenceCriteria(maxIter, confidence))
        
        if result.fitness > fitness:
          fitness = result.fitness
          best_result = result
        count += 1

    return best_result


  def refine_registration(self, source, target, voxel_size, result_ransac, ICPThreshold_factor):
    from open3d import pipelines
    from open3d import utility
    utility.random.seed(99)
    registration = pipelines.registration
    distance_threshold = voxel_size * ICPThreshold_factor
    print(":: Point-to-plane ICP registration is applied on original point")
    print("   clouds to refine the alignment. This time we use a strict")
    print("   distance threshold %.3f." % distance_threshold)
    result = registration.registration_icp(
        source, target, distance_threshold, result_ransac.transformation,
        registration.TransformationEstimationPointToPlane())
    return result


  def takeScreenshot(self,name,description,type=-1):
    # show the message even if not taking a screen shot
    slicer.util.delayDisplay('Take screenshot: '+description+'.\nResult is available in the Annotations module.', 3000)

    lm = slicer.app.layoutManager()
    # switch on the type to get the requested window
    widget = 0
    if type == slicer.qMRMLScreenShotDialog.FullLayout:
      # full layout
      widget = lm.viewport()
    elif type == slicer.qMRMLScreenShotDialog.ThreeD:
      # just the 3D window
      widget = lm.threeDWidget(0).threeDView()
    elif type == slicer.qMRMLScreenShotDialog.Red:
      # red slice window
      widget = lm.sliceWidget("Red")
    elif type == slicer.qMRMLScreenShotDialog.Yellow:
      # yellow slice window
      widget = lm.sliceWidget("Yellow")
    elif type == slicer.qMRMLScreenShotDialog.Green:
      # green slice window
      widget = lm.sliceWidget("Green")
    else:
      # default to using the full window
      widget = slicer.util.mainWindow()
      # reset the type so that the node is set correctly
      type = slicer.qMRMLScreenShotDialog.FullLayout

    # grab and convert to vtk image data
    qimage = ctk.ctkWidgetsUtils.grabWidget(widget)
    imageData = vtk.vtkImageData()
    slicer.qMRMLUtils().qImageToVtkImageData(qimage,imageData)

    annotationLogic = slicer.modules.annotations.logic()
    annotationLogic.CreateSnapShot(name, description, type, 1, imageData)

  def process(self, inputVolume, outputVolume, imageThreshold, invert=False, showResult=True):
    """
    Run the processing algorithm.
    """

    if not inputVolume or not outputVolume:
      raise ValueError("Input or output volume is invalid")

    import time
    startTime = time.time()
    logging.info('Processing started')

    # Compute the thresholded output volume using the "Threshold Scalar Volume" CLI module
    cliParams = {
      'InputVolume': inputVolume.GetID(),
      'OutputVolume': outputVolume.GetID(),
      'ThresholdValue' : imageThreshold,
      'ThresholdType' : 'Above' if invert else 'Below'
      }
    cliNode = slicer.cli.run(slicer.modules.thresholdscalarvolume, None, cliParams, wait_for_completion=True, update_display=showResult)
    # We don't need the CLI module node anymore, remove it to not clutter the scene with it
    slicer.mrmlScene.RemoveNode(cliNode)

    stopTime = time.time()
    logging.info(f'Processing completed in {stopTime-startTime:.2f} seconds')



class BiteIDTest(ScriptedLoadableModuleTest):
  """
    This is the test case for your scripted module.
    Uses ScriptedLoadableModuleTest base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
      """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
      """
    self.setUp()
    self.test_BiteID1()

  def test_BiteID1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")

    # Get/create input data

    import SampleData
    registerSampleData()
    inputVolume = SampleData.downloadSample('TemplateKey1')
    self.delayDisplay('Loaded test data set')

    inputScalarRange = inputVolume.GetImageData().GetScalarRange()
    self.assertEqual(inputScalarRange[0], 0)
    self.assertEqual(inputScalarRange[1], 695)

    outputVolume = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode")
    threshold = 100

    # Test the module logic

    logic = BiteIDLogic()

    # Test algorithm with non-inverted threshold
    logic.process(inputVolume, outputVolume, threshold, True)
    outputScalarRange = outputVolume.GetImageData().GetScalarRange()
    self.assertEqual(outputScalarRange[0], inputScalarRange[0])
    self.assertEqual(outputScalarRange[1], threshold)

    # Test algorithm with inverted threshold
    logic.process(inputVolume, outputVolume, threshold, False)
    outputScalarRange = outputVolume.GetImageData().GetScalarRange()
    self.assertEqual(outputScalarRange[0], inputScalarRange[0])
    self.assertEqual(outputScalarRange[1], inputScalarRange[1])

    self.delayDisplay('Test passed')

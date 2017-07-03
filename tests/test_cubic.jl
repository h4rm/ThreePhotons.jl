#Cubic
doubleRotatedDensity = rotateStructure(rotateStructure(densityCube, 1.5, 1.5, 1.5), -1.5, -1.5, -1.5)
@test similarity(densityCube, doubleRotatedDensity) > 0.85

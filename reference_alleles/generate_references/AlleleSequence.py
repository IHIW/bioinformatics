def featureNames():
    return ['5UTR','Exon1','Intron1','Exon2','Intron2','Exon3','Intron3','Exon4','Intron4','Exon5','Intron5','Exon6','Intron6','Exon7','Intron7','Exon8','3UTR']

class AlleleSequence:
    def __init__(self):
        self.alleleName = None
        self.description = None
        self.featureSequences = {}
        self.accessionNumber = None

    def getLocus(self):
        try:
            currentLocus, nomenclatureFields = self.alleleName.split('*')
            return currentLocus
        except Exception as e:
            return None

    def isFullLength(self):
        featureNames = ','.join(sorted(list(self.featureSequences.keys())))
        #print('These are the feature names:' + str(featureNames))
        # Here's a hack to see if it's full length. I'm sure this is not perfect
        # Only works with 4-8 exons
        return (featureNames in [
            '3UTR,5UTR,Exon1,Exon2,Exon3,Exon4,Intron1,Intron2,Intron3'
            ,'3UTR,5UTR,Exon1,Exon2,Exon3,Exon4,Exon5,Intron1,Intron2,Intron3,Intron4'
            ,'3UTR,5UTR,Exon1,Exon2,Exon3,Exon4,Exon5,Exon6,Intron1,Intron2,Intron3,Intron4,Intron5'
            ,'3UTR,5UTR,Exon1,Exon2,Exon3,Exon4,Exon5,Exon6,Exon7,Intron1,Intron2,Intron3,Intron4,Intron5,Intron6'
            ,'3UTR,5UTR,Exon1,Exon2,Exon3,Exon4,Exon5,Exon6,Exon7,Exon8,Intron1,Intron2,Intron3,Intron4,Intron5,Intron6,Intron7'
            ])

    def getSequence(self):
        completeSequence = ''
        for featureName in featureNames():
            if(featureName in self.featureSequences.keys()):
                completeSequence += self.featureSequences[featureName]
        return completeSequence

from pcraster import *
from pcraster.framework import *
import random
import csv
import numpy


class ModifiedForestFire(DynamicModel):
    def __init__(self):
        DynamicModel.__init__(self)
        setclone('bool200.map')

    def initial(self):
        self.blankMap = readmap("bool200")
        self.fires = readmap("trees200")
        self.tempfires = self.fires
        self.cluster_size_list = []
        self.run = 0
        # Option set for clump function
        setglobaloption("nondiagonal")

    def dynamic(self):
        # Random location for lightning
        random_x = random.randint(0, 199)
        random_y = random.randint(-199, 0)

        # Geometric mass function
        q = 0.001

        f1 = q * (1 - q) ** (1 - 1)
        f2 = q * (1 - q) ** (2 - 1)
        f3 = q * (1 - q) ** (3 - 1)

        f_1 = f1 / (f1 + f2 + f3)
        f_2 = f2 / (f1 + f2 + f3)
        f_3 = f3 / (f1 + f2 + f3)

        neighbourhood = numpy.random.choice([1, 2, 3], p=[f_1, f_2, f_3])
        # <<< Use this for static neighbourhood >>>     neighbourhood = 1
        neighbourhood = int(neighbourhood)

        lightning = ifthenelse(
            pcrand(pcrand(xcoordinate(self.blankMap) > random_x, xcoordinate(self.blankMap) < random_x + 1),
                   pcrand(ycoordinate(self.blankMap) > random_y - 1, ycoordinate(self.blankMap) < random_y)),
            boolean(1), boolean(0))

        # Check wether lightning hit a tree
        hit_tree = ifthenelse(pcrand(lightning, self.fires == 0), boolean(1), boolean(0))
        if maptotal(scalar(hit_tree)) == 1:
            # Use clump function if neighbourhood is 1
            if neighbourhood == 1:
                bool_trees = self.fires == 0
                clump_trees = clump(bool_trees)
                new_burn_cluster = ifthenelse(lightning, clump_trees, nominal(0))
                nbc_total = maptotal(scalar(new_burn_cluster))
                new_burn_cluster = ifthenelse(scalar(clump_trees) == nbc_total, boolean(1), boolean(0))
                cluster_size = maptotal(scalar(new_burn_cluster))
                self.cluster_size_list.append([self.run, int(cluster_size)])
                self.fires = ifthenelse(new_burn_cluster, nominal(1), self.fires)

            else:
                # Density based clustering for neighbourhood 2 or 3
                lightmap = spread(lightning, 0, 1)
                new_burn = ifthenelse(pcror(pcrand(lightmap <= neighbourhood, self.fires == 0), lightning), boolean(1),
                                      boolean(0))
                not_done = True
                while not_done:
                    old_burn = maptotal(scalar(new_burn))
                    old_spread = spread(new_burn, 0, 1)
                    new_burn = ifthenelse(pcrand(old_spread <= neighbourhood, self.fires == 0), boolean(1), boolean(0))
                    if maptotal(scalar(new_burn)) == 0:
                        not_done = False
                    if maptotal(scalar(new_burn)) == old_burn:
                        not_done = False

                # Calculate cluster size
                new_burn = ifthenelse(pcrand(new_burn, self.fires == 0), boolean(1), boolean(0))
                cluster_size = maptotal(scalar(new_burn))
                # Calculate new map
                self.cluster_size_list.append([self.run, int(cluster_size)])
                self.fires = ifthenelse(pcrand(new_burn, self.fires == 0), 1, self.fires)

        # Grow new tree
        tree_growth = 1 / 500
        newTrees = uniform(1) < tree_growth

        self.fires = ifthenelse(pcrand(newTrees, self.fires == 1), nominal(0), self.fires)
        self.run += 1
        print(self.run)

        # Write results to file
        if nrOfTimeSteps <= self.run:
            csvfile = "mffmv001.csv"
            with open(csvfile, "w") as output:
                writer = csv.writer(output)
                writer.writerows(self.cluster_size_list)


nrOfTimeSteps = 50000
myModel = ModifiedForestFire()
dynamicModel = DynamicFramework(myModel, nrOfTimeSteps)
dynamicModel.run()

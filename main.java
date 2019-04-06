import edu.princeton.cs.algs4.StdDraw;
import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.StdRandom;

public class Main {

    public static void main(String[] args) {
        int nCents = 9; //the number of centroids generated
        double[][] cents = randCents(nCents);

        StdDraw.setPenRadius(0.01);
        StdDraw.setPenColor(StdDraw.RED);
        drawCents(cents);

        //setup and call grid-based lloyd's algorithm for centroidal voronoi tesselation problem in 2d between 0 and 1 for xy bounds
        int nIter = 10000;
        int nGrid = 100; //define grid resolution (nGrid x nGrid)

        cents = gridLloyd(cents, nIter, nGrid);

        StdDraw.setPenColor(StdDraw.BLUE);
        drawCents(cents);
        StdDraw.show();
    }

    public static double[][] gridLloyd(double[][] cents, int nIter, int nGrid){
        int nCents = cents.length;
        for (int iter = 0; iter < nIter; iter++) {
            //obtain the distances and assign samples to cluster
            int member[][] = new int[1000][1000];
            double[][] newCents = new double[nCents][2]; //accumulate sum of coordinates in each cluster
            int[] clusterCount = new int[nCents]; //count the number of grid points in each cluster
            double[][] grid = new double[nGrid][nGrid];
            for (int m = 0; m < nGrid; m++) {
                for (int n = 0; n < nGrid; n++) {
                    double d = Integer.MAX_VALUE; //init to large number
                    for (int i = 0; i < nCents; i++) {
                        double d_swap = distance(cents[i][0], cents[i][1], (double) (m) / (double) (nGrid), (double) (n) / (double) (nGrid));
                        if (d_swap < d) { //assign correct distance
                            d = d_swap;
                            member[m][n] = i; //select membership of point for book keeping and visualization :D
                            newCents[i][0] += (double) (m) / (double) (nGrid); //accumulate the sum of the cluster coords
                            newCents[i][1] += (double) (n) / (double) (nGrid);
                            clusterCount[i]++; //keep count of points in each cluster
                        }
                    }
                }
            }
            for (int i = 0; i < nCents; i++) { //average points within clusters to assign new centroids
                cents[i][0] = newCents[i][0] / clusterCount[i];
                cents[i][1] = newCents[i][1] / clusterCount[i];
                StdOut.printf("Iter #%d: Centroid %d, (%f, %f)\n", iter, i, cents[i][0], cents[i][1]);
            }
        }
        return cents;
    }
    //generate and return N uniformly distributed points between [0,1)
    public static double[][] randCents(int N) {
        double[][] cents;
        cents = new double[N][2];
        for (int i = 0; i < N; i++) {
            cents[i][0] = StdRandom.uniform();
            cents[i][1] = StdRandom.uniform();
        }
        return cents;
    }

    //Call StdDraw.point on each centroid (Does not call StdDraw.show(), however)
    public static void drawCents(double[][] cents) {
        for (int i = 0; i < cents.length; i++) {
            StdDraw.point(cents[i][0], cents[i][1]);
        }
    }

    //Euclidean Distance Function
    public static double distance(double x1, double y1, double x2, double y2) {
        return Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2));
    }
}

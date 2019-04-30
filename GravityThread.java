
/*
 * "Physics" part of code adapted from Dan Schroeder's applet at:
 *
 *     http://physics.weber.edu/schroeder/software/mdapplet.html
 */
import java.awt.*;
import java.util.concurrent.BrokenBarrierException;
import java.util.concurrent.CyclicBarrier;
import javax.swing.*;

public class GravityThread extends Thread {

    final static int THREADS = 8;
    static CyclicBarrier Barrier = new CyclicBarrier(THREADS);
    static Display display = new Display();
    // Size of simulation

    final static int N = 2000;  // Number of "stars"
    final static double BOX_WIDTH = 100.0;

    // Initial state
    final static double RADIUS = 20.0;  // of randomly populated sphere

    final static double ANGULAR_VELOCITY = 0.4;
    // controls total angular momentum

    // Simulation
    final static double DT = 0.002;  // Time step

    // Display
    final static int WINDOW_SIZE = 800;
    final static int DELAY = 0;
    final static int OUTPUT_FREQ = 5;

    // Star positions
    volatile static double[] x = new double[N];
    volatile static double[] y = new double[N];
    volatile static double[] z = new double[N];

    // Star velocities
    volatile static double[] vx = new double[N];
    volatile static double[] vy = new double[N];
    volatile static double[] vz = new double[N];

    // Star accelerations
    volatile static double[] ax = new double[N];
    volatile static double[] ay = new double[N];
    volatile static double[] az = new double[N];

    public static long startTime; // will be called by threads later on

    public static void main(String args[]) throws Exception {

        // generate stars in the galaxy initial pos and angular velocity 
    	GenerateGalaxy();

        // create threads
        GravityThread[] threads = new GravityThread[THREADS];
        for (int i = 0; i < THREADS; i++) {
            threads[i] = new GravityThread(i);
            threads[i].start();
        }

        startTime = System.currentTimeMillis();

        for (int i = 0; i < THREADS; i++) {
            threads[i].join();
        }

        long elapsedTime = System.currentTimeMillis() - startTime;
        System.out.println("Multithread Elapsed Time = " + elapsedTime + "ms");
        System.exit(0); // lazy me forgets to close the jframe
    }

    int me;

    public GravityThread(int me) {
        this.me = me;
    }

    public void Sleep() {
        try {
            Thread.sleep(DELAY);
        } catch (InterruptedException ex) {
           System.exit(1);
        }
    }
    
    static public void Sync() {
        try {
            Barrier.await();
        } catch (BrokenBarrierException | InterruptedException ex) {
            System.exit(1);
        }
    }
    
    @Override
    public void run() {
        int iter = 0;
        // the N is going to be a new value,
        // we need to divide it by number of threads
        int b = N / THREADS;
        int Begin = me * b ;
        int End = Begin + b ;
        
        // try sync after position updates as different stages require surrounding parts
        System.out.println("thread " + me + " start: " + Begin + " end: " + End);
        
        while (iter < 2000) { 
            if (me == 0 && iter % OUTPUT_FREQ == 0) {
                //System.out.println("iter = " + iter + ", time = " + iter * DT);
                display.repaint();
            }
            
            // Verlet integration:
            // http://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
            double dtOver2 = 0.5 * DT;
            double dtSquaredOver2 = 0.5 * DT * DT;
            for (int i = Begin; i < End; i++) {
                // update position
                x[i] += (vx[i] * DT) + (ax[i] * dtSquaredOver2);
                y[i] += (vy[i] * DT) + (ay[i] * dtSquaredOver2);
                z[i] += (vz[i] * DT) + (az[i] * dtSquaredOver2);
                // update velocity halfway
                vx[i] += (ax[i] * dtOver2);
                vy[i] += (ay[i] * dtOver2);
                vz[i] += (az[i] * dtOver2);
            }
            Sync();
            
            // need to look at the axis in this as well
            computeAccelerations(Begin, End); 
            Sync();
            
            // update velocities with new acceleration
            for (int i = Begin; i < End; i++) {
                vx[i] += (ax[i] * dtOver2);
                vy[i] += (ay[i] * dtOver2);
                vz[i] += (az[i] * dtOver2);
            }

            Sync();
            
            iter++;
           // System.out.println("thread " + me + " iter = " + iter);
            
        }
        
        //System.out.println("thread " + me + " Elapsed Time = " + (System.currentTimeMillis() - startTime) + "ms");
    }

    // Compute accelerations of all stars from current positions:
    static void computeAccelerations(int begin, int end) {

        double dx, dy, dz;  // separations in x and y directions
        double dx2, dy2, dz2, rSquared, r, rCubedInv, fx, fy, fz;

        // no point reseting everything
        for (int i = begin; i < end; i++) {
            ax[i] = 0.0;
            ay[i] = 0.0;
            az[i] = 0.0;
        }

        // Interaction forces (gravity)
        for (int i = begin; i < end; i++) {
            for (int j = 0; j < N; j++) {  // loop over all distinct pairs
                    
            	if (i != j) // Required  or sets location to 0
                    break;
                // Vector version of inverse square law
                dx = x[i] - x[j];
                dy = y[i] - y[j];
                dz = z[i] - z[j];
                dx2 = dx * dx;
                dy2 = dy * dy;
                dz2 = dz * dz;
                rSquared = dx2 + dy2 + dz2;
                r = Math.sqrt(rSquared);
                rCubedInv = 1.0 / (rSquared * r);
                fx = -rCubedInv * dx;
                fy = -rCubedInv * dy;
                fz = -rCubedInv * dz;

                ax[i] += fx;  // add this force on to i's acceleration (mass = 1)
                ay[i] += fy;
                az[i] += fz;
               // ax[j] -= fx;  // Newton's 3rd law
              //  ay[j] -= fy;
              //  az[j] -= fz;
            }  
        }
    }

    // Intitial generation of galaxy 
    static void GenerateGalaxy() {
    	 // Define initial state of stars
        // Randomly choose plane for net angular velocity
        //nx, ny, nz represents angular vel (or rotation speed)
        double nx = 2 * Math.random() - 1;
        double ny = 2 * Math.random() - 1;
        double nz = 2 * Math.random() - 1;
        double norm = 1.0 / Math.sqrt(nx * nx + ny * ny + nz * nz);
        nx *= norm;
        ny *= norm;
        nz *= norm;

        for (int i = 0; i < N; i++) {

            // Place star randomly in sphere of specified radius
            double rx, ry, rz, r;
            do {
                rx = (2 * Math.random() - 1) * RADIUS;
                ry = (2 * Math.random() - 1) * RADIUS;
                rz = (2 * Math.random() - 1) * RADIUS;
                r = Math.sqrt(rx * rx + ry * ry + rz * rz);
            } while (r > RADIUS);

            x[i] = 0.5 * BOX_WIDTH + rx;
            y[i] = 0.5 * BOX_WIDTH + ry;
            z[i] = 0.5 * BOX_WIDTH + rz;

            vx[i] = ANGULAR_VELOCITY * (ny * rz - nz * ry);
            vy[i] = ANGULAR_VELOCITY * (nz * rx - nx * rz);
            vz[i] = ANGULAR_VELOCITY * (nx * ry - ny * rx);
        }
    }
    
    static class Display extends JPanel {

        static final double SCALE = WINDOW_SIZE / BOX_WIDTH;

        Display() {

            setPreferredSize(new Dimension(WINDOW_SIZE, WINDOW_SIZE));

            JFrame frame = new JFrame("MD");
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setContentPane(this);
            frame.pack();
            frame.setVisible(true);
        }

        public void paintComponent(Graphics g) {
            g.setColor(Color.BLACK);
            g.fillRect(0, 0, WINDOW_SIZE, WINDOW_SIZE);
            g.setColor(Color.WHITE);
            for (int i = 0; i < N; i++) {
            	
            	// choose colour based on thread
            	g.setColor(Color.WHITE);
            	//g.setColor(GetThreadColour(i));
            	               
                int gx = (int) (SCALE * x[i]);
                int gy = (int) (SCALE * y[i]);
                if (0 <= gx && gx < WINDOW_SIZE && 0 < gy && gy < WINDOW_SIZE) {
                    g.fillRect(gx, gy, 1, 1);
                }
            }
        }
        
        private Color GetThreadColour(int i) {
        	double threadNo = i / (N / (double)THREADS);
        	if (threadNo <= 1)
                return Color.RED;
        	else if (threadNo <= 2)
        		return Color.WHITE;
        	else if (threadNo <= 3)
        		return Color.GREEN;
        	else if (threadNo <= 4)
        		return Color.BLUE;
        	else if (threadNo <= 5)
        		return Color.ORANGE;
        	else if (threadNo <= 6)
        		return Color.PINK;
        	else if (threadNo <= 7)
        		return Color.MAGENTA;
        	else if (threadNo <= 8)
        		return Color.CYAN;
        	else
                throw new UnsupportedOperationException();
        }
    }
}

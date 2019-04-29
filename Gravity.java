
/*
 * "Physics" part of code adapted from Dan Schroeder's applet at:
 *
 *     http://physics.weber.edu/schroeder/software/mdapplet.html
 */

import java.awt.* ;
import javax.swing.* ;

public class Gravity {

    // Size of simulation

    final static int N = 2000 ;  // Number of "stars"
    final static double BOX_WIDTH = 100.0 ;


    // Initial state

    final static double RADIUS = 20.0 ;  // of randomly populated sphere

    final static double ANGULAR_VELOCITY = 0.4 ;
           // controls total angular momentum


    // Simulation

    final static double DT = 0.002 ;  // Time step


    // Display

    final static int WINDOW_SIZE = 1000 ;
    final static int DELAY = 0 ;
    final static int OUTPUT_FREQ = 5 ;


    // Star positions
    static double [] x = new double [N] ;
    static double [] y = new double [N] ;
    static double [] z = new double [N] ;

    // Star velocities
    static double [] vx = new double [N] ;
    static double [] vy = new double [N] ;
    static double [] vz = new double [N] ;

    // Star accelerations
    static double [] ax = new double [N] ;
    static double [] ay = new double [N] ;
    static double [] az = new double [N] ;


    public static void main(String args []) throws Exception {

        Display display = new Display() ;

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

        

        // ... or just rotate in x, y plane
        //double nx = 0, ny = 0, nz = 1.0 ;

        // ... or just rotate in x, z plane
        //double nx = 0, ny = 1.0, nz = 0 ;

        for(int i = 0 ; i < N ; i++) {

            // Place star randomly in sphere of specified radius
            double rx, ry, rz, r ;
            do {
                rx = (2 * Math.random() - 1) * RADIUS ;
                ry = (2 * Math.random() - 1) * RADIUS ;
                rz = (2 * Math.random() - 1) * RADIUS ;
                r = Math.sqrt(rx * rx + ry * ry + rz * rz) ;
            } while(r > RADIUS) ;

            x [i] = 0.5 * BOX_WIDTH + rx ;
            y [i] = 0.5 * BOX_WIDTH + ry ;
            z [i] = 0.5 * BOX_WIDTH + rz ;

            vx [i] = ANGULAR_VELOCITY * (ny * rz - nz * ry) ; 
            vy [i] = ANGULAR_VELOCITY * (nz * rx - nx * rz) ; 
            vz [i] = ANGULAR_VELOCITY * (nx * ry - ny * rx) ; 
        }
         

        int iter = 0 ;
        long startTime = System.currentTimeMillis();
        while(iter < 2000) {

            if(iter % OUTPUT_FREQ == 0) {
                System.out.println("iter = " + iter + ", time = " + iter * DT) ;
                display.repaint() ;
                Thread.sleep(DELAY) ;
            }

            // Verlet integration:
            // http://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet

            double dtOver2 = 0.5 * DT;
            double dtSquaredOver2 = 0.5 * DT * DT;
            for (int i = 0; i < N; i++) {
                // update position
                x[i] += (vx[i] * DT) + (ax[i] * dtSquaredOver2);
                y[i] += (vy[i] * DT) + (ay[i] * dtSquaredOver2);
                z[i] += (vz[i] * DT) + (az[i] * dtSquaredOver2);
                // update velocity halfway
                vx[i] += (ax[i] * dtOver2);
                vy[i] += (ay[i] * dtOver2);
                vz[i] += (az[i] * dtOver2);
            }

            computeAccelerations();

            for (int i = 0; i < N; i++) {
                // finish updating velocity with new acceleration
                vx[i] += (ax[i] * dtOver2);
                vy[i] += (ay[i] * dtOver2);
                vz[i] += (az[i] * dtOver2);
            }

            iter++ ;
        }
        
        long elapsedTime = System.currentTimeMillis() - startTime;
        System.out.println("ElapsedTime = " + elapsedTime + "ms");
    }

    // Compute accelerations of all stars from current positions:
    static void computeAccelerations() {

        double dx, dy, dz;  // separations in x and y directions
        double dx2, dy2, dz2, rSquared, r, rCubedInv, fx, fy, fz;

        for (int i = 0; i < N; i++) {
            ax[i] = 0.0;
            ay[i] = 0.0;
            az[i] = 0.0;
        }

        // Interaction forces (gravity)
        // This is where the program spends most of its time.

        // (NOTE: use of Newton's 3rd law below to essentially half number
        // of calculations needs some care in a parallel version.
        // A naive decomposition on the i loop can lead to a race condition
        // because you are assigning to ax[j], etc.
        // You can remove these assignments and extend the j loop to a fixed
        // upper bound of N, or, for extra credit, find a cleverer solution!)

        for (int i = 1; i < N; i++) {
            for (int j = 0; j < i; j++) {  // loop over all distinct pairs

                // Vector version of inverse square law
                dx = x[i] - x[j];
                dy = y[i] - y[j];
                dz = z[i] - z[j];
                dx2 = dx * dx;
                dy2 = dy * dy;
                dz2 = dz * dz;
                rSquared = dx2 + dy2 + dz2 ;
                r = Math.sqrt(rSquared) ;
                rCubedInv = 1.0 / (rSquared * r) ;
                fx = - rCubedInv * dx;
                fy = - rCubedInv * dy;
                fz = - rCubedInv * dz;

                ax[i] += fx;  // add this force on to i's acceleration (mass = 1)
                ay[i] += fy;
                az[i] += fz;
                ax[j] -= fx;  // Newton's 3rd law
                ay[j] -= fy;
                az[j] -= fz;
            }
        }
    }


    static class Display extends JPanel {

        static final double SCALE = WINDOW_SIZE / BOX_WIDTH ;

        Display() {

            setPreferredSize(new Dimension(WINDOW_SIZE, WINDOW_SIZE)) ;

            JFrame frame = new JFrame("MD");
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setContentPane(this);
            frame.pack();
            frame.setVisible(true);
        }

        public void paintComponent(Graphics g) {
            g.setColor(Color.BLACK) ;
            g.fillRect(0, 0, WINDOW_SIZE, WINDOW_SIZE) ;
            g.setColor(Color.WHITE) ;
            for(int i = 0 ; i < N ; i++) {
                int gx = (int) (SCALE * x [i]) ;
                int gy = (int) (SCALE * y [i]) ;
                if(0 <= gx && gx < WINDOW_SIZE && 0 < gy && gy < WINDOW_SIZE) { 
                    g.fillRect(gx, gy, 1, 1) ;
                }
            } 
        }
    }
}
    
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gravityaparapi;

import com.aparapi.Kernel;
import com.aparapi.Range;
import com.aparapi.device.Device;

import java.awt.*;
import java.awt.Color;
import java.util.concurrent.CyclicBarrier;
import javax.swing.*;

/**
 *
 * @author up788458
 */
public class GravityAparapi {

    final static int THREADS = 8;
    static CyclicBarrier Barrier = new CyclicBarrier(THREADS);
    // Size of simulation

    final static int N = 6000;  // Number of "stars"
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
    final static int OUTPUT_FREQ = 25;

    // non final static variables are not supported, will have to init and make during execution
    // Star positions
//    static double[] x = new double[N];
//    static double[] y = new double[N];
//    static double[] z = new double[N];
//
//    // Star velocities
//    static double[] vx = new double[N];
//    static double[] vy = new double[N];
//    static double[] vz = new double[N];
//
//    // Star accelerations
//    static double[] ax = new double[N];
//    static double[] ay = new double[N];
//    static double[] az = new double[N];
    public static void main(String args[]) throws Exception {
        
        // Star positions
        double[] x = new double[N];
        double[] y = new double[N];
        double[] z = new double[N];

        // Star velocities
        double[] vx = new double[N];
        double[] vy = new double[N];
        double[] vz = new double[N];

        // Star accelerations
        double[] ax = new double[N];
        double[] ay = new double[N];
        double[] az = new double[N];
        
        // generate stars in the galaxy initial pos and angular velocity
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
        
        Display display = new Display();
        
        Device gpu = Device.best();
        System.out.println(gpu.toString());

        // create kernel here
        Kernel kernel = new Kernel() {
            @Override
            public void run() {
                int globalid = getGlobalId();

                double dx, dy, dz;  // separations in x and y directions
                double dx2, dy2, dz2, rSquared, r, rCubedInv, fx, fy, fz;

                // Interaction forces (gravity)
                for (int j = 0; j < N; j++) {  // loop over all distinct pairs

                    // does not support break >:(
                    //if (globalid == j) 
                    //    break;
                    if (globalid != j) {
                        // Vector version of inverse square law
                        dx = x[globalid] - x[j];
                        dy = y[globalid] - y[j];
                        dz = z[globalid] - z[j];
                        dx2 = dx * dx;
                        dy2 = dy * dy;
                        dz2 = dz * dz;
                        rSquared = dx2 + dy2 + dz2;
                        r = Math.sqrt(rSquared);
                        rCubedInv = 1.0 / (rSquared * r);
                        fx = -rCubedInv * dx;
                        fy = -rCubedInv * dy;
                        fz = -rCubedInv * dz;

                        ax[globalid] = fx;  // add this force on to i's acceleration (mass = 1)
                        ay[globalid] = fy;
                        az[globalid] = fz;
                        // ax[j] -= fx;  // Newton's 3rd law
                        //  ay[j] -= fy;
                        //  az[j] -= fz;
                    }
                }

            }
        ;

        };                
        
        //Range range = Range.create(N);
        Range range = gpu.createRange(N);

        long startTime = System.currentTimeMillis();
               // Thread.sleep(10000);
        
        int iter = 0;
        while (iter < 2000) {
            if (iter % OUTPUT_FREQ == 0) {
                System.out.println("iter = " + iter + ", time = " + iter * DT);
                display.UpdatePos(x, y);
                display.repaint();
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

            // need to look at the axis in this as well
            //computeAccelerations();
            // a kernel for each i iteration
            // gid represents what i is
            // replaces compute acceleratinos mthod
            kernel.execute(range);

            // update velocities with new acceleration
            for (int i = 0; i < N; i++) {
                vx[i] += (ax[i] * dtOver2);
                vy[i] += (ay[i] * dtOver2);
                vz[i] += (az[i] * dtOver2);
            }

            iter++;
        }

        long elapsedTime = System.currentTimeMillis() - startTime;

        System.out.println("GPU Elapsed Time = " + elapsedTime + "ms");
        System.exit(0); // lazy me forgets to close the jframe
    }

    static class Display extends JPanel {

        static final double SCALE = WINDOW_SIZE / BOX_WIDTH;
        double[] x = new double[N];
        double[] y = new double[N];
        
        Display() {
            for (int i = 0; i < N; i++){             
                x[i] = 0;
                y[i] = 0;
            }

            setPreferredSize(new Dimension(WINDOW_SIZE, WINDOW_SIZE));

            JFrame frame = new JFrame("MD");
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setContentPane(this);
            frame.pack();
            frame.setVisible(true);
        }
        
        public void UpdatePos(double[] x, double[] y){
            this.x = x;
            this.y = y;
        }

        @Override
        public void paintComponent(Graphics g) {
            g.setColor(Color.BLACK);
            g.fillRect(0, 0, WINDOW_SIZE, WINDOW_SIZE);
            g.setColor(Color.WHITE);
            for (int i = 0; i < N; i++) {

                // choose colour based on thread
                // comment out to keep white
                //g.setColor(GetThreadColour(i));
                
                // code cannot be fetched from non static field, will nedd to update      
              
                    try {
                        int gx = (int) (SCALE * x[i]);
                        int gy = (int) (SCALE * y[i]);

                        if (0 <= gx && gx < WINDOW_SIZE && 0 < gy && gy < WINDOW_SIZE) 
                            g.fillRect(gx, gy, 1, 1);
                    }
                    catch (NullPointerException e) {
                        System.out.println(i + " is null! WTF");
                        continue;
                    }
                    
                
                
            }
        }

        private Color GetThreadColour(int i) {
            double threadNo = i / (N / (double) THREADS);
            if (threadNo <= 1) {
                return Color.RED;
            } else if (threadNo <= 2) {
                return Color.WHITE;
            } else if (threadNo <= 3) {
                return Color.GREEN;
            } else if (threadNo <= 4) {
                return Color.BLUE;
            } else if (threadNo <= 5) {
                return Color.ORANGE;
            } else if (threadNo <= 6) {
                return Color.PINK;
            } else if (threadNo <= 7) {
                return Color.MAGENTA;
            } else if (threadNo <= 8) {
                return Color.CYAN;
            } else {
                throw new UnsupportedOperationException();
            }
        }
    }

}

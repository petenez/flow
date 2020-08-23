// Author: Petri Hirvonen, petenez@gmail.com, 2 September 2019

package flowplotter;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import static java.lang.System.exit;

public class FlowPlotter {

	public static void main(String[] args) throws IOException {

		String input1;
		String input2;
		String output;
		int Nx;	// width of the system in pixels
		int Ny;	// height
		// parameters for optional looping over time steps
		int start;
		int incr;
		int end;
		if(args.length != 3 && args.length != 6) {
			System.out.println("Error: Invalid syntax.");
			exit(1);
		}
		if(args.length == 6) {
			start = Integer.parseInt(args[3]);
			incr = Integer.parseInt(args[4]);
			end = Integer.parseInt(args[5]);
			if(start > end || incr < 1) {
				System.out.println("Error: Invalid batch mode parameters.");
				exit(1);
			}
		} else {
			start = 0;
			incr = 1;
			end = 0;
		}

		// loop over different time steps
		for(int m = start; m <= end; m += incr) {
			input1 = args[0];
			input2 = args[1];
			output = args[2];
			if(args.length == 6) {	// loop over different time steps
				String[] pieces = input1.trim().split("#");	// split 'input' at "#"
				input1 = "";	// input emptied
				if(pieces.length >= 1) input1 += pieces[0];	// append piece before "#"
				input1 += m;	// append time step
				if(pieces.length >= 2) input1 += pieces[1];	// append piece after "#"

				pieces = input2.trim().split("#");	// split 'input' at "#"
				input2 = "";	// input emptied
				if(pieces.length >= 1) input2 += pieces[0];	// append piece before "#"
				input2 += m;	// append time step
				if(pieces.length >= 2) input2 += pieces[1];	// append piece after "#"

				pieces = output.trim().split("#");	// repeat for 'output'
				output = "";
				if(pieces.length >= 1) output += pieces[0];
				output += m;
				if(pieces.length >= 2) output += pieces[1];
			}

			File f = new File(input1);
			if(!f.exists()) {	// check if input file is found
				if(start < end) {
					System.out.println("Warning: Input file " + input1 + " not found. Continuing iteration.");
					continue;
				}
				else {
					System.out.println("Error: Input file " + input1 + " not found.");
					exit(1);
				}
			}
			f = new File(input2);
			if(!f.exists()) {	// check if input file is found
				if(start < end) {
					System.out.println("Warning: Input file " + input2 + " not found. Continuing iteration.");
					continue;
				}
				else {
					System.out.println("Error: Input file " + input2 + " not found.");
					exit(1);
				}
			}

			BufferedReader reader1 = new BufferedReader(new FileReader(input1));
			BufferedReader reader2 = new BufferedReader(new FileReader(input2));
			String[] temp = reader1.readLine().trim().split("\\s+");
			reader2.readLine();
			Nx = Integer.parseInt(temp[0]);
			Ny = Integer.parseInt(temp[1]);
			double[][] ndata = new double[Ny][Nx];	// data array
			double[][] vxdata = new double[Ny][Nx];	// data array
			double[][] vydata = new double[Ny][Nx];	// data array

			// determines extrema
			double nmin = 1.0e300;	// initially huge just to be safe
			double nmax = -1.0e300;	// initially -huge ...
			double vmax = 0.0;

			double n, vx, vy, v;
			for(int j = 0; j < Ny; j++) {
				for(int i = 0; i < Nx; i++) {
					n = Double.parseDouble(reader1.readLine());
					ndata[j][i] = n;
					if(n < nmin) nmin = n;	// check if current minimum
					if(n > nmax) nmax = n;	// maximum
					temp = reader2.readLine().trim().split("\\s+");
					vx = Double.parseDouble(temp[0]);
					vy = Double.parseDouble(temp[1]);
					vxdata[j][i] = vx;
					vydata[j][i] = vy;
					v = Math.sqrt(vx*vx + vy*vy);
					if(v > vmax) vmax = v;
				}
			}

			reader1.close();
			reader2.close();

			// sets pixels and writes output
			BufferedImage image = new BufferedImage(Nx, Ny, BufferedImage.TYPE_INT_RGB);
			double _vmax = 1.0/vmax;
			double _2pi = 0.5/Math.PI;
			for(int j = 0; j < Ny; j++) {
				for(int i = 0; i < Nx; i++) {
					float b = (float)((ndata[j][i] - nmin)/(nmax -nmin));
					float s = (float)(Math.sqrt(vxdata[j][i]*vxdata[j][i] + vydata[j][i]*vydata[j][i])*_vmax);
					float h = (float)(_2pi*Math.atan2(vydata[j][i], vxdata[j][i]) + 0.5);
					image.setRGB(i, Ny - j - 1, Color.HSBtoRGB(h, s, b));
				}
			}
			ImageIO.write(image, "png", new File(output));
		}
	}
}
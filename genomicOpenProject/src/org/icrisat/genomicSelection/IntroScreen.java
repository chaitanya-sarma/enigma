/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.icrisat.genomicSelection;

import java.awt.SplashScreen;
import java.awt.Graphics2D;
import java.awt.Font;
import java.awt.Dimension;
import java.awt.Color;
import java.awt.geom.Rectangle2D;

/**
 *
 * @author Mohan
 */
//SplashScreen class used for the startup of the application 
public class IntroScreen { //class to display loading libraies requried

    static SplashScreen mySplash;
    static Graphics2D splashGraphics;
    static Rectangle2D.Double splashTextArea;
    static Rectangle2D.Double splashProgressArea;
    static Font font;

    //#####################################################################
    //main method of the IntroScreen class
    public static void main(String[] args) {
        splashInit();
        appInit();
    }
    //##########################################################################
    //method to show updating text on splashScreen

    private static void appInit() {
        for (int i = 1; i <= 10; i++) {
            int pctDone = i * 10;
            splashText("Loading modules...please wait");
            splashProgress(pctDone);
        }
        try {
            Thread.sleep(1000);
        } catch (InterruptedException ie) {
            System.out.println(ie);
        }
     
        //Initilization of Main class GsMain
        GsMain gselection = new GsMain();
        gselection.startup();
    }

    //#########################################################################
    //method to set the demissions of the splash screen image  
    private static void splashInit() {
        mySplash = SplashScreen.getSplashScreen();
        if (mySplash != null) {
            Dimension ssDim = mySplash.getSize();
            int height = ssDim.height;
            int width = ssDim.width;
            splashTextArea = new Rectangle2D.Double(15., height * 0.88, width * .45, 32.);
            splashProgressArea = new Rectangle2D.Double(width * .55, height * .92, width * .4, 12);
            splashGraphics = mySplash.createGraphics();
            font = new Font("Dialog", Font.PLAIN, 14);
            splashGraphics.setFont(font);
            splashText("Starting");
            splashProgress(0);

        }
    }

    //#################################################################################
    //setting the text to on the splash screen
    public static void splashText(String str) {
        if (mySplash != null && mySplash.isVisible()) {
            splashGraphics.setPaint(Color.LIGHT_GRAY);
            splashGraphics.fill(splashTextArea);
            splashGraphics.setPaint(Color.BLACK);
            splashGraphics.drawString(str, (int) (splashTextArea.getX() + 10), (int) (splashTextArea.getY() + 15));
            mySplash.update();
        }
    }
//########################################################################################
    //setting the progress number 

    public static void splashProgress(int pct) {
        if (mySplash != null && mySplash.isVisible()) {
            splashGraphics.setPaint(Color.BLUE);
            splashGraphics.fill(splashTextArea);
            splashGraphics.setPaint(Color.WHITE);
            splashGraphics.fill(splashProgressArea);

            int x = (int) splashProgressArea.getMinX();
            int y = (int) splashProgressArea.getMinY();
            int wid = (int) splashProgressArea.getWidth();
            int hgt = (int) splashProgressArea.getHeight();
            int doneWidth = Math.round(pct * wid / 100.f);
            doneWidth = Math.max(0, Math.min(doneWidth, wid - 1));
            splashGraphics.setPaint(Color.blue);
            splashGraphics.fillRect(x, y + 1, doneWidth, hgt - 1);
        }
    }
}

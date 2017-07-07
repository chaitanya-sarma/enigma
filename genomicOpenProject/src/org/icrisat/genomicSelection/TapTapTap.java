package org.icrisat.genomicSelection;

import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.AlphaComposite;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JLayer;
import javax.swing.Timer;
import javax.swing.ButtonGroup;
import javax.swing.JRadioButton;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JTextField;
import javax.swing.JComponent;
import javax.swing.plaf.LayerUI;


/*
 * class which is used in showing the progress display on dialog box of 
 * Data summary and progressDialog class 
 */
//start region 
//not used region 
public class TapTapTap {

    public static void main(String[] args) {
        javax.swing.SwingUtilities.invokeLater(new Runnable() {

            @Override
            public void run() {
                new TapTapTap().createUI();
            }
        });
    }
    private JButton mOrderButton;

    public void createUI() {
        JFrame f = new JFrame("TapTapTap");

        final WaitLayerUI layerUI = new WaitLayerUI();
        JPanel panel = createPanel();
        JLayer<JPanel> jlayer = new JLayer<>(panel, layerUI);

        final Timer stopper = new Timer(40000, new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent ae) {
                layerUI.stop();
            }
        });
        stopper.setRepeats(false);

        mOrderButton.addActionListener(
                new ActionListener() {

                    @Override
                    public void actionPerformed(ActionEvent ae) {
                        layerUI.start();
                        if (!stopper.isRunning()) {
                            stopper.start();
                        }
                    }
                });

        f.add(jlayer);

        f.setSize(300, 200);
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.setLocationRelativeTo(null);
        f.setVisible(true);
    }

    private JPanel createPanel() {
        JPanel p = new JPanel();

        ButtonGroup entreeGroup = new ButtonGroup();
        JRadioButton radioButton;
        p.add(radioButton = new JRadioButton("Beef", true));
        entreeGroup.add(radioButton);
        p.add(radioButton = new JRadioButton("Chicken"));
        entreeGroup.add(radioButton);
        p.add(radioButton = new JRadioButton("Vegetable"));
        entreeGroup.add(radioButton);

        p.add(new JCheckBox("Ketchup"));
        p.add(new JCheckBox("Mustard"));
        p.add(new JCheckBox("Pickles"));

        p.add(new JLabel("Special requests:"));
        p.add(new JTextField(20));

        mOrderButton = new JButton("Place Order");
        p.add(mOrderButton);

        return p;
    }
}
//end region 

//start region
/*
 * class used for displaying progress in DataSummary class
 */
class WaitLayerUI extends LayerUI<JPanel> implements ActionListener {

    private boolean mIsRunning;
    private boolean mIsFadingOut;
    private Timer mTimer;
    private int mAngle;
    private int mFadeCount;
    private int mFadeLimit = 15;

    @Override
    public void paint(Graphics g, JComponent c) {
        int w=345 ;  //for setting the x axis in display in data summary
        int h=100;   //for setting y axis in display in datasummary
        // Paint the view.
        super.paint(g, c);

        if (!mIsRunning) {
            return;
        }

        Graphics2D g2 = (Graphics2D) g.create();

        float fade = (float) mFadeCount / (float) mFadeLimit;
        // Gray it out.
        // Paint the wait indicator.
        int s = Math.min(110,110) / 5;
        int cx = w / 1;
        int cy = h / 1;
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
        g2.setStroke(new BasicStroke(s / 4, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
        g2.setPaint(Color.blue);
        g2.rotate(Math.PI * mAngle / 180, cx, cy);
        for (int i = 0; i < 12; i++) {
            float scale = (11.0f - (float) i) / 11.0f;
            g2.drawLine(cx + s, cy, cx + s * 2, cy);
            g2.rotate(-Math.PI / 6, cx, cy);
            g2.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, scale * fade));
        }
        g2.dispose();
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        if (mIsRunning) {
            firePropertyChange("tick", 0, 1);
            mAngle += 3;
            if (mAngle >= 360) {
                mAngle = 0;
            }
            if (mIsFadingOut) {
                if (--mFadeCount == 0) {
                    mIsRunning = false;
                    mTimer.stop();
                }
            } else if (mFadeCount < mFadeLimit) {
                mFadeCount++;
            }
        }
    }

    public void start() {
        if (mIsRunning) {
            return;
        }

        // Run a thread for animation.
        mIsRunning = true;
        mIsFadingOut = false;
        mFadeCount = 0;
        int fps = 24;
        int tick = 1000 / fps;
        mTimer = new Timer(tick, this);
        mTimer.start();
    }

    public void stop() {
        mIsFadingOut = true;
    }

    @Override
    public void applyPropertyChange(PropertyChangeEvent pce, JLayer l) {
        if ("tick".equals(pce.getPropertyName())) {
            l.repaint();
        }
    }
}
//end region 

//start region 
/*
 * class used for progressDisplay class 
*/
class WaitLayerUI1 extends LayerUI<JPanel> implements ActionListener {

    private boolean mIsRunning;
    private boolean mIsFadingOut;
    private Timer mTimer;
    private int mAngle;
    private int mFadeCount;
    private int mFadeLimit = 15;

    @Override
    public void paint(Graphics g, JComponent c) {
        int w=185;
        int h=160;
        // Paint the view.
        super.paint(g, c);

        if (!mIsRunning) {
            return;
        }

        Graphics2D g2 = (Graphics2D) g.create();

        float fade = (float) mFadeCount / (float) mFadeLimit;
      // Paint the wait indicator.
        int s = Math.min(80,80) / 5;
        int cx = w / 1;
        int cy = h / 1;
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
        g2.setStroke(new BasicStroke(s / 4, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
        g2.setPaint(Color.blue);
        g2.rotate(Math.PI * mAngle / 180, cx, cy);
        for (int i = 0; i < 12; i++) {
            float scale = (11.0f - (float) i) / 11.0f;
            g2.drawLine(cx + s, cy, cx + s * 2, cy);
            g2.rotate(-Math.PI / 6, cx, cy);
            g2.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, scale * fade));
        }
        g2.dispose();
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        if (mIsRunning) {
            firePropertyChange("tick", 0, 1);
            mAngle += 3;
            if (mAngle >= 360) {
                mAngle = 0;
            }
            if (mIsFadingOut) {
                if (--mFadeCount == 0) {
                    mIsRunning = false;
                    mTimer.stop();
                }
            } else if (mFadeCount < mFadeLimit) {
                mFadeCount++;
            }
        }
    }

    public void start() {
        if (mIsRunning) {
            return;
        }
        // Run a thread for animation.
        mIsRunning = true;
        mIsFadingOut = false;
        mFadeCount = 0;
        int fps = 24;
        int tick = 1000 / fps;
        mTimer = new Timer(tick, this);
        mTimer.start();
    }

    public void stop() {
        mIsFadingOut = true;
    }

    @Override
    public void applyPropertyChange(PropertyChangeEvent pce, JLayer l) {
        if ("tick".equals(pce.getPropertyName())) {
            l.repaint();
        }
    }
}
//end region 

package comp559.lcp;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;

import javax.vecmath.GMatrix;
import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;

/**
 * Implementation of a contact constraint.
 * @author kry
 */
public class Contact {

    /** Next available contact index, used for determining which rows of the jacobian a contact uses */
    static public int nextContactIndex = 0;
    
    /** Index of this contact, determines its rows in the jacobian */
    int index;
    
    /** First RigidBody in contact */
    RigidBody body1;
    
    /** Second RigidBody in contact */
    RigidBody body2;
    
    /** Contact normal in world coordinates */
    Vector2d normal = new Vector2d();
    
    /** Position of contact point in world coordinates */
    Point2d contactW = new Point2d();


    double[] jacobianRowOne;
    double[] jacobianRowTwo;
    double[] massMat;
    
    /**
     * Creates a new contact, and assigns it an index
     * @param body1
     * @param body2
     * @param contactW
     * @param normal
     */
    public Contact( RigidBody body1, RigidBody body2, Point2d contactW, Vector2d normal ) {
        this.body1 = body1;
        this.body2 = body2;
        this.contactW.set( contactW );
        this.normal.set( normal );        
        index = nextContactIndex++;        
        // TODO: Objective 3 - Compute and store the contact Jacobian
        Vector2d body1ToContact = new Vector2d(contactW.x - body1.x.x, contactW.y - body1.x.y);
        Vector2d body2ToContact = new Vector2d( contactW.x - body2.x.x, contactW.y - body2.x.y);

        Vector2d tangent = new Vector2d(-normal.y, normal.x);
        jacobianRowOne = new double[] { -normal.x, -normal.y, body1ToContact.y * normal.x - body1ToContact.x * normal.y,
                                            normal.x, normal.y, body2ToContact.x * normal.y - body2ToContact.y * normal.x};
        jacobianRowTwo = new double[] { -tangent.x, -tangent.y, body1ToContact.y * tangent.x - body1ToContact.x * tangent.y,
                                            tangent.x, tangent.y, body2ToContact.x * tangent.y - body2ToContact.y * tangent.x};
        massMat = new double[] {body1.minv, body1.minv, body1.jinv, body2.minv, body2.minv, body2.jinv};
    }
    
    /**
     * Draws the contact points
     * @param drawable
     */
    public void display( GLAutoDrawable drawable ) {
        GL2 gl = drawable.getGL().getGL2();
        gl.glPointSize(3);
        gl.glColor3f(.7f,0,0);
        gl.glBegin( GL.GL_POINTS );
        gl.glVertex2d(contactW.x, contactW.y);
        gl.glEnd();
    }
    
    /**
     * Draws the connections between bodies to visualize the 
     * the adjacency structure of the matrix as a graph.
     * @param drawable
     */
    public void displayConnection( GLAutoDrawable drawable ) {
        GL2 gl = drawable.getGL().getGL2();
        // draw a line between the two bodies but only if they're both not pinned
        if ( !body1.pinned && ! body2.pinned ) {
            gl.glLineWidth(2);
            gl.glColor4f(0,.3f,0, 0.5f);
            gl.glBegin( GL.GL_LINES );
            gl.glVertex2d(body1.x.x, body1.x.y);
            gl.glVertex2d(body2.x.x, body2.x.y);
            gl.glEnd();
        }
    }
    
}

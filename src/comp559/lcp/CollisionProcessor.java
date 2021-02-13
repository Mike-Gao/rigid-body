package comp559.lcp;

import java.util.*;

import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;

import mintools.parameters.BooleanParameter;
import mintools.parameters.DoubleParameter;
import mintools.parameters.IntParameter;
import mintools.swing.CollapsiblePanel;
import mintools.swing.VerticalFlowPanel;

/**
 * Class for detecting and resolving collisions.  Currently this class uses penalty forces between rigid bodies.
 * @author kry
 */
public class CollisionProcessor {

    private List<RigidBody> bodies;
    
    /**
     * The current contacts that resulted in the last call to process collisions
     */
    public ArrayList<Contact> contacts = new ArrayList<Contact>();
    
    /**
     * Creates this collision processor with the provided set of bodies
     * @param bodies
     */
    public CollisionProcessor( List<RigidBody> bodies ) {
        this.bodies = bodies;
    }
    
    /** keeps track of the time used for collision detection on the last call */
    double collisionDetectTime = 0;
    
    /** keeps track of the time used to solve the LCP based velocity update on the last call */
    double collisionSolveTime = 0;

    public Map<Contact, double[]> pair = new HashMap<>();
    
    /**
     * Processes all collisions 
     * @param dt time step
     */
    public void processCollisions( double dt ) {
        contacts.clear();
        Contact.nextContactIndex = 0;
        
        long now = System.nanoTime();
        broadPhase();
        collisionDetectTime = ( System.nanoTime() - now ) * 1e-9;
                
        if ( contacts.size() > 0  && doLCP.getValue() ) {
            now = System.nanoTime();
            double mu = friction.getValue();

            // TODO: Objective 3 - Compute velocity update with iterative solve of contact constraint matrix.
            int bodyNum = bodies.size();
            int contactNum = contacts.size();

            // construct d, the diagonal of A matrix and the b matrix
            double[] dMat = new double[2 * contactNum];
            double[] bMat = new double[2 * contactNum];
            for (Contact cnct : contacts) {
                // for dmat
                double a = 0;
                double b = 0;
                // ---------------------------------------- //
                // for bmat
                double normalConstraint = 0;
                double frictionConstraint = 0;
                double restitutionParam = 0;
                double[] u = new double[] {cnct.body1.v.x, cnct.body1.v.y, cnct.body1.omega,
                        cnct.body2.v.x, cnct.body2.v.y, cnct.body2.omega};

                double[] f = new double[] {cnct.body1.force.x, cnct.body1.force.y, cnct.body1.torque,
                        cnct.body2.force.x, cnct.body2.force.y, cnct.body2.torque};
                // ---------------------------------------- //
                for (int i = 0; i < 6; i++) {
                    // for dMat
                    a += cnct.jacobianRowOne[i] * cnct.jacobianRowOne[i] * cnct.massMat[i];
                    b += cnct.jacobianRowTwo[i] * cnct.jacobianRowTwo[i] * cnct.massMat[i];
                    // for bMat
                    restitutionParam += cnct.jacobianRowOne[i] * u[i] * restitution.getValue();
                    normalConstraint += cnct.jacobianRowOne[i] * (u[i] + dt * f[i] * cnct.massMat[i]);
                    frictionConstraint += cnct.jacobianRowTwo[i] * (u[i] + dt * f[i] * cnct.massMat[i]);
                }
                dMat[2 * cnct.index + 0] = a;
                dMat[2 * cnct.index + 1] = b;
                bMat[2 * cnct.index + 0] = normalConstraint + restitutionParam;
                bMat[2 * cnct.index + 1] = frictionConstraint;

            }

            double[] bPrimeMat = new double[bMat.length];
            for (int i = 0; i < bPrimeMat.length; i++) {
                bPrimeMat[i] = bMat[i] / dMat[i];
            }

            double[] lambdaI = new double[2 * contactNum];
            double[] friction = new double[contactNum];
            double[] deltaV = new double[3 * bodyNum];

            // construct d and b matrix

            for (int iter = 0; iter < iterations.getValue(); iter++) {
                // TODO: Objective 4 - Optimization
                if (useShuffle.getValue()) Collections.shuffle(contacts);
                for (Contact cnct : contacts) {
                    // TODO: Objective 5 - Warm Start
                    if (warmStart.getValue()) {
                        if(pair.containsKey(cnct)) {
                            lambdaI[2*cnct.index] = pair.get(cnct)[0];
                            lambdaI[2*cnct.index] = pair.get(cnct)[1];
                        }
                    }
                    double lambdaNormal = lambdaI[2 * cnct.index] - bPrimeMat[2 * cnct.index];

                    for(int i = 0; i < 3; i++) {
                        lambdaNormal -= cnct.jacobianRowOne[i + 0] / dMat[2 * cnct.index] * deltaV[3 * cnct.body1.index + i];
                        lambdaNormal -= cnct.jacobianRowOne[i + 3] / dMat[2 * cnct.index] * deltaV[3 * cnct.body2.index + i];
                    }

                    // make sure lambdaNormal > 0
                    lambdaNormal = Math.max(0, lambdaNormal);

                    friction[cnct.index] = mu * lambdaNormal;

                    double delta = lambdaNormal - lambdaI[2 * cnct.index];

                    lambdaI[2 * cnct.index] = lambdaNormal;

                    double[] colIOfT = new double[6];
                    for (int i = 0; i < 6; i++){
                        colIOfT[i] = cnct.massMat[i] * cnct.jacobianRowOne[i];
                    }

                    for (int i = 0; i < 3; i++) {
                        deltaV[3 * cnct.body1.index + i] += colIOfT[i] * delta;
                        deltaV[3 * cnct.body2.index + i] += colIOfT[i + 3] * delta;
                    }

                    double lambdaFriction = lambdaI[2 * cnct.index + 1] - bPrimeMat[2 * cnct.index + 1];
                    for (int i = 0; i < 3; i++) {
                        lambdaFriction -= deltaV[3 * cnct.body1.index + i] * cnct.jacobianRowTwo[i] / dMat[2 * cnct.index + 1];
                        lambdaFriction -= deltaV[3 * cnct.body2.index + i] * cnct.jacobianRowTwo[i + 3] / dMat[2 * cnct.index + 1];
                    }

                    lambdaFriction = Math.max(lambdaFriction, -friction[cnct.index]);
                    lambdaFriction = Math.min(lambdaFriction, friction[cnct.index]);

                    delta = lambdaFriction - lambdaI[2 * cnct.index + 1];
                    lambdaI[2 * cnct.index + 1] = lambdaFriction;

                    for (int i = 0; i < colIOfT.length; i++) {
                        colIOfT[i] = cnct.massMat[i] * cnct.jacobianRowTwo[i];
                    }

                    for (int i = 0; i < 3; i++) {
                        deltaV[3 * cnct.body1.index + i] += colIOfT[i] * delta;
                        deltaV[3 * cnct.body2.index + i] += colIOfT[i + 3] * delta;
                    }

                }
            }

            for (RigidBody rb : bodies) {
                rb.v.x += deltaV[3 * rb.index];
                rb.v.y += deltaV[3 * rb.index + 1];
                rb.omega += deltaV[3 * rb.index + 2];
            }
            if (warmStart.getValue()){
                for(Contact cnct : contacts){
                    pair.put(cnct, new double[] {lambdaI[2* cnct.index], lambdaI[2* cnct.index +1]});
                }
            }

            collisionSolveTime = (System.nanoTime() - now) * 1e-9;
        }
    }
    
    /**
     * Checks for collisions between bodies.  Note that you can optionaly implement some broad
     * phase test such as spatial hashing to reduce the n squared body-body tests.
     * Currently this does the naive n squared collision check.
     */
    private void broadPhase() {
        // Naive n squared body test.. might not be that bad for small number of bodies 
        visitID++;
        for ( RigidBody b1 : bodies ) {
            for ( RigidBody b2 : bodies ) { // not so inefficient given the continue on the next line
                if ( b1.index >= b2.index ) continue;
                if ( b1.pinned && b2.pinned ) continue;                
                narrowPhase( b1, b2 );                
            }
        }        
    }
    
    /**
     * Checks for collision between boundary blocks on two rigid bodies.
     * TODO: Objective 2 - This needs to be improved as the n-squared block test is too slow!
     * @param body1
     * @param body2
     */
    private void narrowPhase( RigidBody body1, RigidBody body2 ) {
        if ( ! useBVTree.getValue() ) {
            for ( Block b1 : body1.blocks ) {
                for ( Block b2 : body2.blocks ) {
                    processCollision( body1, b1, body2, b2 );
                }
            }
        } else {
            // TODO: Objective 2 - Implement code to use hierarchical collision detection on body pairs
            BVNode bvn1 = body1.root;
            BVNode bvn2 = body2.root;
            bvn1.boundingDisc.updatecW();
            bvn2.boundingDisc.updatecW();

            if (bvn1.boundingDisc.intersects(bvn2.boundingDisc)) {
                // we traverse the tree here
                traverseBVH(body1, bvn1, body2, bvn2);
            }
        }
    }

    private void updateVisitID(BVNode in, int id){
        if (in.visitID != id) {
            in.boundingDisc.updatecW();
            in.visitID = id;
        }
    }

    private void traverseBVH(RigidBody rb1, BVNode bvn1, RigidBody rb2, BVNode bvn2){
        if (bvn1.isLeaf() && bvn2.isLeaf()) {
            processCollision(rb1, bvn1.leafBlock, rb2, bvn2.leafBlock);
        } else if (bvn1.isLeaf() && !bvn2.isLeaf()) {
            BVNode left = bvn2.child1;
            BVNode right = bvn2.child2;
            updateVisitID(left, this.visitID);
            updateVisitID(right, this.visitID);
            if (bvn1.boundingDisc.intersects(left.boundingDisc)) {
                traverseBVH(rb1, bvn1, rb2, left);
            }
            if (bvn1.boundingDisc.intersects(right.boundingDisc)) {
                traverseBVH(rb1, bvn1, rb2, right);
            }
        } else if (!bvn1.isLeaf() && bvn2.isLeaf()) {
            BVNode left = bvn1.child1;
            BVNode right = bvn1.child2;
            updateVisitID(left, this.visitID);
            updateVisitID(right, this.visitID);
            if (bvn2.boundingDisc.intersects(left.boundingDisc)) {
                traverseBVH(rb1, left, rb2, bvn2);
            }
            if (bvn2.boundingDisc.intersects(right.boundingDisc)) {
                traverseBVH(rb1, right, rb2, bvn2);
            }
        } else {
            BVNode left, right;
            if (bvn1.boundingDisc.r > bvn2.boundingDisc.r) {
                left = bvn1.child1;
                right = bvn1.child2;
                updateVisitID(left, this.visitID);
                updateVisitID(right, this.visitID);
                if (bvn2.boundingDisc.intersects(left.boundingDisc)) {
                    traverseBVH(rb1, left, rb2, bvn2);
                }
                if (bvn2.boundingDisc.intersects(right.boundingDisc)) {
                    traverseBVH(rb1, right, rb2, bvn2);
                }
            } else {
                left = bvn2.child1;
                right = bvn2.child2;
                updateVisitID(left, this.visitID);
                updateVisitID(right, this.visitID);
                if (bvn1.boundingDisc.intersects(left.boundingDisc)) {
                    traverseBVH(rb1, bvn1, rb2, left);
                }
                if (bvn1.boundingDisc.intersects(right.boundingDisc)) {
                    traverseBVH(rb1, bvn1, rb2, right);
                }
            }
        }

    }


    /** 
     * The visitID is used to tag boundary volumes that are visited in 
     * a given time step.  Marking boundary volume nodes as visited during
     * a time step allows for a visualization of those used, but it can also
     * be used to more efficiently update the centeres of bounding volumes
     * (i.e., call a BVNode's updatecW method at most once on any given timestep)
     */
    int visitID = 0;
    
    /**
     * Resets the state of the collision processor by clearing all
     * currently identified contacts, and reseting the visitID for
     * tracking the bounding volumes used
     */
    public void reset() {
        contacts.clear();
        Contact.nextContactIndex = 0;
        visitID = 0;            
    }
    
    // some working variables for processing collisions
    private Point2d tmp1 = new Point2d();
    private Point2d tmp2 = new Point2d();
    private Point2d contactW = new Point2d();
    private Vector2d force = new Vector2d();
    private Vector2d contactV1 = new Vector2d();
    private Vector2d contactV2 = new Vector2d();
    private Vector2d relativeVelocity = new Vector2d();
    private Vector2d normal = new Vector2d();
        
    /**
     * Processes a collision between two bodies for two given blocks that are colliding.
     * Currently this implements a penalty force
     * @param body1
     * @param b1
     * @param body2
     * @param b2
     */
    private void processCollision( RigidBody body1, Block b1, RigidBody body2, Block b2 ) {        
        double k = contactSpringStiffness.getValue();
        double c1 = contactSpringDamping.getValue();
        double threshold = separationVelocityThreshold.getValue();
        boolean useSpring = enableContactSpring.getValue();
        boolean useDamping = enableContactDamping.getValue();
        
        body1.transformB2W.transform( b1.pB, tmp1 );
        body2.transformB2W.transform( b2.pB, tmp2 );
        double distance = tmp1.distance(tmp2);
        if ( distance < Block.radius * 2 ) {
            // contact point at halfway between points 
            // NOTE: this assumes that the two blocks have the same radius!
            contactW.interpolate( tmp1, tmp2, .5 );
            // contact normal
            normal.sub( tmp2, tmp1 );
            normal.normalize();
            // create the contact
            Contact contact = new Contact( body1, body2, contactW, normal);
            // simple option... add to contact list...
            contacts.add( contact );
            if ( ! doLCP.getValue() ) {
                // compute relative body velocity at contact point
                body1.getSpatialVelocity( contactW, contactV1 );
                body2.getSpatialVelocity( contactW, contactV2 );
                relativeVelocity.sub( contactV1, contactV2 );
                if ( -relativeVelocity.dot( normal ) < threshold ) {
                    if ( useSpring ) {
                        // spring force
                        double interpenetration = distance - Block.radius * 2; // a negative quantity
                        force.scale( -interpenetration * k, normal );
                        body2.applyContactForceW(contactW, force);
                        force.scale(-1);
                        body1.applyContactForceW(contactW, force);
                    }
                    if ( useDamping ) {
                        // spring damping forces!
                        // vertical
                        force.scale( relativeVelocity.dot(normal) * c1, normal );                    
                        body2.applyContactForceW( contactW, force );
                        force.scale(-1);
                        body1.applyContactForceW( contactW, force );
                    }
                }
            }
        }
    }
   
    /** Stiffness of the contact penalty spring */
    private DoubleParameter contactSpringStiffness = new DoubleParameter("penalty contact stiffness", 1e3, 1, 1e5 );
    
    /** Viscous damping coefficient for the contact penalty spring */
    private DoubleParameter contactSpringDamping = new DoubleParameter("penalty contact damping", 10, 1, 1e4 );
    
    /** Threshold for the relative velocity in the normal direction, for determining if spring force will be applied. */
    private DoubleParameter separationVelocityThreshold = new DoubleParameter( "penalty separation velocity threshold (controls bounce)", 1e-9, 1e-9, 1e3 );
    
    /** Enables the contact penalty spring */
    private BooleanParameter enableContactSpring = new BooleanParameter("enable penalty contact spring", true );
    
    /** Enables damping of the contact penalty spring */
    private BooleanParameter enableContactDamping = new BooleanParameter("enable penalty contact damping", true );
    
    /** Restitution parameter for contact constraints */
    public DoubleParameter restitution = new DoubleParameter( "restitution (bounce)", 0, 0, 1 );
    
    /** Coulomb friction coefficient for contact constraint */
    public DoubleParameter friction = new DoubleParameter("Coulomb friction", 0.33, 0, 2 );
    
    /** Number of iterations to use in projected Gauss Seidel solve */
    public IntParameter iterations = new IntParameter("iterations for GS solve", 10, 1, 500);
    
    /** Flag for switching between penalty based contact and contact constraints */
    private BooleanParameter doLCP = new BooleanParameter( "do LCP solve", false );
    
    /** Flag for enabling the use of hierarchical collision detection for body pairs */
    private BooleanParameter useBVTree = new BooleanParameter( "use BVTree", false );

    /** Flag for enabling randomization */
    private BooleanParameter useShuffle = new BooleanParameter("use Shuffle", false);

    private BooleanParameter warmStart = new BooleanParameter("use Warm Start", false);

    /**
     * @return controls for the collision processor
     */
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
        vfp.setBorder( new TitledBorder("Collision Processing Controls") );
        vfp.add( useBVTree.getControls() );
        vfp.add( useShuffle.getControls());
        vfp.add( doLCP.getControls() );
        vfp.add( warmStart.getControls());
        vfp.add( iterations.getSliderControls() );
        vfp.add( restitution.getSliderControls(false) );
        vfp.add( friction.getSliderControls(false) );
        
        VerticalFlowPanel vfp2 = new VerticalFlowPanel();
        vfp2.setBorder( new TitledBorder("penalty method controls") );
        vfp2.add( contactSpringStiffness.getSliderControls(true) );
        vfp2.add( contactSpringDamping.getSliderControls(true) );
        vfp2.add( separationVelocityThreshold.getSliderControls( true ) );
        vfp2.add( enableContactDamping.getControls() );
        vfp2.add( enableContactSpring.getControls() );
        
        CollapsiblePanel cp = new CollapsiblePanel(vfp2.getPanel());
        cp.collapse();
        vfp.add( cp );        
        return vfp.getPanel();
    }
    
}

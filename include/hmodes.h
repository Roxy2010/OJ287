/*
 * l=2 harmonics of h
 */

static void eval_R(REAL t, const REAL* f, const CBwaveODE& ode,
                         const ObserverParameters& obs, REAL* R)
{
    // relative position vector
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);
    REAL n[3] = {rx/r, ry/r, rz/r};
    REAL dm = ode.m1 - ode.m2;

    // relative velocity vector
    REAL v[3] = {eval_vx(t, f, ode, obs),
                 eval_vy(t, f, ode, obs),
                 eval_vz(t, f, ode, obs)};
    REAL vsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    REAL rdot = (rx*v[0] + ry*v[1] + rz*v[2])/r;

    REAL eta = ode.eta;
    REAL Gm_r = ode.m/r;

    // spins
    REAL cSs1 = ode.m1*ode.m1;
    REAL S1x = cSs1*f[I_s1x];
    REAL S1y = cSs1*f[I_s1y];
    REAL S1z = cSs1*f[I_s1z];
    REAL cSs2 = ode.m2*ode.m2;
    REAL S2x = cSs2*f[I_s2x];
    REAL S2y = cSs2*f[I_s2y];
    REAL S2z = cSs2*f[I_s2z];
    REAL Sx = S1x + S2x;
    REAL Sy = S1y + S2y;
    REAL Sz = S1z + S2z;

    REAL Deltax = ode.m*(S2x/ode.m2 - S1x/ode.m1);
    REAL Deltay = ode.m*(S2y/ode.m2 - S1y/ode.m1);
    REAL Deltaz = ode.m*(S2z/ode.m2 - S1z/ode.m1);
    REAL S[3], Delta[3];
    S[0] = Sx;
    S[1] = Sy;
    S[2] = Sz;
    Delta[0] = Deltax;
    Delta[1] = Deltay;
    Delta[2] = Deltaz;
    REAL n_cross_v[3];
    n_cross_v[0] = n[1]*v[2] - n[2]*v[1];
    n_cross_v[1] = n[2]*v[0] - n[0]*v[2];
    n_cross_v[2] = n[0]*v[1] - n[1]*v[0];
    REAL P15QSO_A = 0;
    for(int i = 0; i < 3; ++i) {
        P15QSO_A += n_cross_v[i]*(12*S[i] + 6*dm/ode.m*Delta[i]);
    }
    REAL v_cross_B[3];
    v_cross_B[0] = v[1]*(9*Sz + 5*dm/ode.m*Deltaz) - v[2]*(9*Sy + 5*dm/ode.m*Deltay);
    v_cross_B[1] = v[2]*(9*Sx + 5*dm/ode.m*Deltax) - v[0]*(9*Sz + 5*dm/ode.m*Deltaz);
    v_cross_B[2] = v[0]*(9*Sy + 5*dm/ode.m*Deltay) - v[1]*(9*Sx + 5*dm/ode.m*Deltax);
    REAL n_cross_D[3];
    n_cross_D[0] = n[1]*(2*Sz + 2*dm/ode.m*Deltaz) - n[2]*(2*Sy + 2*dm/ode.m*Deltay);
    n_cross_D[1] = n[2]*(2*Sx + 2*dm/ode.m*Deltax) - n[0]*(2*Sz + 2*dm/ode.m*Deltaz);
    n_cross_D[2] = n[0]*(2*Sy + 2*dm/ode.m*Deltay) - n[1]*(2*Sx + 2*dm/ode.m*Deltax);
    REAL n_cross_E[3];
    n_cross_E[0] = n[1]*(12*Sz + 6*dm/ode.m*Deltaz) - n[2]*(12*Sy + 6*dm/ode.m*Deltay);
    n_cross_E[1] = n[2]*(12*Sx + 6*dm/ode.m*Deltax) - n[0]*(12*Sz + 6*dm/ode.m*Deltaz);
    n_cross_E[2] = n[0]*(12*Sy + 6*dm/ode.m*Deltay) - n[1]*(12*Sx + 6*dm/ode.m*Deltax);

    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
            REAL sum = 0;
            if((ode.hterms & H_Q) != 0) {
                REAL Q = 2*(v[i]*v[j] - Gm_r*n[i]*n[j]);
                sum += Q;
            }
            if((ode.hterms & H_PQ) != 0) {
                REAL PQ = 2/3.0*ode.m/r*rdot*(5 + 3*eta)*(n[i]*v[j] + n[j]*v[i])
                    + ((1-3*eta)*vsq - 2/3.0*(2 - 3*eta)*ode.m/r)*v[i]*v[j]
                    + ode.m/r*((1-3*eta)*rdot*rdot - 1/3.0*(10 + 3*eta)*vsq
                            + 29/3.0*ode.m/r)*n[i]*n[j] ;
                sum += PQ;
            }
            if((ode.hterms & H_P15QSO) != 0) {
                REAL P15QSO = 1.0/rsq
                    *( 2*n[i]*n[j]*P15QSO_A - n[i]*v_cross_B[j] - n[j]*v_cross_B[i]
                            - (v[i]*n_cross_D[j] + v[j]*n_cross_D[i])
                            + rdot*(n[i]*n_cross_E[j] + n[j]*n_cross_E[i]) );
                sum += P15QSO;
            }
// Todo: Tail and SS
          R[3*i + j] = 2*ode.mu*sum;
        }
    }
}

static void eval_Pn(REAL t, const REAL* f, const CBwaveODE& ode,
                         const ObserverParameters& obs, REAL* Pn)
{
    // relative position vector
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);
    REAL n[3] = {rx/r, ry/r, rz/r};

    // relative velocity vector
    REAL v[3] = {eval_vx(t, f, ode, obs),
                 eval_vy(t, f, ode, obs),
                 eval_vz(t, f, ode, obs)};

    REAL vsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    REAL rdot = (rx*v[0] + ry*v[1] + rz*v[2])/r;

    REAL eta = ode.eta;
    REAL Gm_r = ode.m/r;
    REAL dm = ode.m1 - ode.m2;

    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
            REAL sum = 0;
            if((ode.hterms & H_P05Q) != 0) {
                REAL P05Q = dm/ode.m*3*Gm_r*(n[i]*v[j] + n[j]*v[i]
                                          - rdot*n[i]*n[j]);
                sum += P05Q;
            }
            if((ode.hterms & H_P15Q) != 0) {
                REAL P15Q = dm/(12.0*r)*( n[i]*n[j]*rdot*(rdot*rdot*(15 - 90*eta)
                                - vsq*(63 - 54*eta)
                                + Gm_r*(242 - 24*eta))
                            - rdot*v[i]*v[j]*(186 + 24*eta)
                            + (n[i]*v[j] + n[j]*v[i])
                            *(rdot*rdot*(63 + 54*eta)
                                - Gm_r*(128 - 36*eta) + vsq*(33 - 18*eta)) );
                sum += P15Q;
            }
          Pn[3*i + j] = 2*ode.mu*sum;
        }
    }
}

static void eval_Pv(REAL t, const REAL* f, const CBwaveODE& ode,
                         const ObserverParameters& obs, REAL* Pv)
{
    // relative position vector
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);
    REAL n[3] = {rx/r, ry/r, rz/r};

    // relative velocity vector
    REAL v[3] = {eval_vx(t, f, ode, obs),
                 eval_vy(t, f, ode, obs),
                 eval_vz(t, f, ode, obs)};

    REAL vsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    REAL rdot = (rx*v[0] + ry*v[1] + rz*v[2])/r;

    REAL eta = ode.eta;
    REAL Gm_r = ode.m/r;
    REAL dm = ode.m1 - ode.m2;

    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
            REAL sum = 0;
            if((ode.hterms & H_P05Q) != 0) {
                REAL P05Q = dm/ode.m*(Gm_r*n[i]*n[j] - 2*v[i]*v[j]);
                sum += P05Q;
            }
            if((ode.hterms & H_P15Q) != 0) {
                REAL P15Q = dm/ode.m*( 0.5*v[i]*v[j]*(Gm_r*(3 - 8*eta) - 2*vsq*(1 - 5*eta))
                            - 0.5*(n[i]*v[j] + n[j]*v[i])*Gm_r*rdot*(7 + 4*eta)
                            - n[i]*n[j]*Gm_r*(3/4.0*(1 - 2*eta)*rdot*rdot
                                + 1/3.0*(26 - 3*eta)*Gm_r
                                - 1/4.0*(7 - 2*eta)*vsq) );
                sum += P15Q;
            }
          Pv[3*i + j] = 2*ode.mu*sum;
        }
    }
}

static void eval_Pnn(REAL t, const REAL* f, const CBwaveODE& ode,
                         const ObserverParameters& obs, REAL* Pnn)
{
    // relative position vector
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);
    REAL n[3] = {rx/r, ry/r, rz/r};

    // relative velocity vector
    REAL v[3] = {eval_vx(t, f, ode, obs),
                 eval_vy(t, f, ode, obs),
                 eval_vz(t, f, ode, obs)};
    REAL vsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    REAL rdot = (rx*v[0] + ry*v[1] + rz*v[2])/r;

    REAL eta = ode.eta;

    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
            REAL sum = 0;
            if((ode.hterms & H_PQ) != 0) {
                REAL PQ = (1-3*eta)/3.0*ode.m/r*((3*vsq - 15*rdot*rdot + 7*ode.m/r)*n[i]*n[j]
                                + 15*rdot*(n[i]*v[j] + n[j]*v[i])
                                - 14*v[i]*v[j]);
                sum += PQ;
            }
          Pnn[3*i + j] = 2*ode.mu*sum;
        }
    }
}

static void eval_Pnv(REAL t, const REAL* f, const CBwaveODE& ode,
                         const ObserverParameters& obs, REAL* Pnv)
{
    // relative position vector
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);
    REAL n[3] = {rx/r, ry/r, rz/r};

    // relative velocity vector
    REAL v[3] = {eval_vx(t, f, ode, obs),
                 eval_vy(t, f, ode, obs),
                 eval_vz(t, f, ode, obs)};
    REAL rdot = (rx*v[0] + ry*v[1] + rz*v[2])/r;

    REAL eta = ode.eta;

    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
            REAL sum = 0;
            if((ode.hterms & H_PQ) != 0) {
                REAL PQ = (1-3*eta)/3.0*4*ode.m/r*(3*rdot*n[i]*n[j]
                        - 4*(n[i]*v[j] + n[j]*v[i]));
                sum += PQ;
            }
          Pnv[3*i + j] = 2*ode.mu*sum;
        }
    }
}

static void eval_Pvv(REAL t, const REAL* f, const CBwaveODE& ode,
                         const ObserverParameters& obs, REAL* Pvv)
{
    // relative position vector
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);
    REAL n[3] = {rx/r, ry/r, rz/r};

    // relative velocity vector
    REAL v[3] = {eval_vx(t, f, ode, obs),
                 eval_vy(t, f, ode, obs),
                 eval_vz(t, f, ode, obs)};

    REAL eta = ode.eta;

    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
            REAL sum = 0;
            if((ode.hterms & H_PQ) != 0) {
                REAL PQ = (1-3*eta)/3.0*2*(3*v[i]*v[j] - ode.m/r*n[i]*n[j]);
                sum += PQ;
            }
          Pvv[3*i + j] = 2*ode.mu*sum;
        }
    }
}

static void eval_Pnnn(REAL t, const REAL* f, const CBwaveODE& ode,
                         const ObserverParameters& obs, REAL* Pnnn)
{
    // relative position vector
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);
    REAL n[3] = {rx/r, ry/r, rz/r};

    // relative velocity vector
    REAL v[3] = {eval_vx(t, f, ode, obs),
                 eval_vy(t, f, ode, obs),
                 eval_vz(t, f, ode, obs)};

    REAL vsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    REAL rdot = (rx*v[0] + ry*v[1] + rz*v[2])/r;

    REAL eta = ode.eta;
    REAL Gm_r = ode.m/r;
    REAL dm = ode.m1 - ode.m2;

    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
            REAL sum = 0;
            if((ode.hterms & H_P15Q) != 0) {
                REAL P15Q = dm*(1-2*eta)/r
                    *( 5/4.0*(3*vsq - 7*rdot*rdot + 6*Gm_r)*rdot*n[i]*n[j]
                         - 1/12.0*(21*vsq - 105*rdot*rdot + 44*Gm_r)
                                *(n[i]*v[j] + n[j]*v[i])
                                - 17/2.0*rdot*v[i]*v[j] );
                sum += P15Q;
            }
          Pnnn[3*i + j] = 2*ode.mu*sum;
        }
    }
}

static void eval_Pnnv(REAL t, const REAL* f, const CBwaveODE& ode,
                         const ObserverParameters& obs, REAL* Pnnv)
{
    // relative position vector
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);
    REAL n[3] = {rx/r, ry/r, rz/r};

    // relative velocity vector
    REAL v[3] = {eval_vx(t, f, ode, obs),
                 eval_vy(t, f, ode, obs),
                 eval_vz(t, f, ode, obs)};

    REAL vsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    REAL rdot = (rx*v[0] + ry*v[1] + rz*v[2])/r;

    REAL eta = ode.eta;
    REAL Gm_r = ode.m/r;
    REAL dm = ode.m1 - ode.m2;

    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
            REAL sum = 0;
            if((ode.hterms & H_P15Q) != 0) {
                REAL P15Q = dm*(1-2*eta)*0.25/r
                    *( (45*rdot*rdot - 9*vsq - 28*Gm_r)*n[i]*n[j]
                                + 58*v[i]*v[j]
                                - 54*rdot*(n[i]*v[j] + n[j]*v[i]) );
                sum += P15Q;
            }
          Pnnv[3*i + j] = 2*ode.mu*sum;
        }
    }
}

static void eval_Pnvv(REAL t, const REAL* f, const CBwaveODE& ode,
                         const ObserverParameters& obs, REAL* Pnvv)
{
    // relative position vector
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);
    REAL n[3] = {rx/r, ry/r, rz/r};

    // relative velocity vector
    REAL v[3] = {eval_vx(t, f, ode, obs),
                 eval_vy(t, f, ode, obs),
                 eval_vz(t, f, ode, obs)};

//    REAL vsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    REAL rdot = (rx*v[0] + ry*v[1] + rz*v[2])/r;

    REAL eta = ode.eta;
//    REAL Gm_r = ode.m/r;
    REAL dm = ode.m1 - ode.m2;

    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
            REAL sum = 0;
            if((ode.hterms & H_P15Q) != 0) {
                REAL P15Q = dm*(1-2*eta)*1.5/r
                    *( 5*(n[i]*v[j] + n[j]*v[i]) - 3*rdot*n[i]*n[j] );
                sum += P15Q;
            }
          Pnvv[3*i + j] = 2*ode.mu*sum;
        }
    }
}

static void eval_Pvvv(REAL t, const REAL* f, const CBwaveODE& ode,
                         const ObserverParameters& obs, REAL* Pvvv)
{
    // relative position vector
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);
    REAL n[3] = {rx/r, ry/r, rz/r};

    // relative velocity vector
    REAL v[3] = {eval_vx(t, f, ode, obs),
                 eval_vy(t, f, ode, obs),
                 eval_vz(t, f, ode, obs)};

//    REAL vsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
//    REAL rdot = (rx*v[0] + ry*v[1] + rz*v[2])/r;

    REAL eta = ode.eta;
    REAL Gm_r = ode.m/r;
    REAL dm = ode.m1 - ode.m2;

    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
            REAL sum = 0;
            if((ode.hterms & H_P15Q) != 0) {
                REAL P15Q = dm/ode.m*(1-2*eta)*0.5
                    *( Gm_r*n[i]*n[j] - 4*v[i]*v[j] );
                sum += P15Q;
            }
          Pvvv[3*i + j] = 2*ode.mu*sum;
        }
    }
}

static FUNCDEF(hp22) {
    REAL Rij[9];
    REAL Pnij[9];
    REAL Pvij[9];
    REAL Pnnij[9];
    REAL Pnvij[9];
    REAL Pvvij[9];
    REAL Pnnnij[9];
    REAL Pnnvij[9];
    REAL Pnvvij[9];
    REAL Pvvvij[9];
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);
    REAL nx = rx/r;  
    REAL ny = ry/r;  
    REAL nz = rz/r;
   
    REAL vx = eval_vx(t, f, ode, obs);
    REAL vy = eval_vy(t, f, ode, obs);    
    REAL vz = eval_vz(t, f, ode, obs);  

    REAL rdot = (rx*vx + ry*vy + rz*vz)/r;
    REAL dm = ode.m1 - ode.m2;

    REAL cSs1 = ode.m1*ode.m1;
    REAL S1x = cSs1*f[I_s1x];
    REAL S1y = cSs1*f[I_s1y];
    REAL S1z = cSs1*f[I_s1z];
    REAL cSs2 = ode.m2*ode.m2;
    REAL S2x = cSs2*f[I_s2x];
    REAL S2y = cSs2*f[I_s2y];
    REAL S2z = cSs2*f[I_s2z];
    REAL Sx = S1x + S2x;
    REAL Sy = S1y + S2y;
    REAL Sz = S1z + S2z;
    REAL Dx = ode.m*(S2x/ode.m2 - S1x/ode.m1);
    REAL Dy = ode.m*(S2y/ode.m2 - S1y/ode.m1);
    REAL Dz = ode.m*(S2z/ode.m2 - S1z/ode.m1);

    eval_R(t, f, ode, obs, Rij);
    eval_Pn(t, f, ode, obs, Pnij);
    eval_Pv(t, f, ode, obs, Pvij);    
    eval_Pnn(t, f, ode, obs, Pnnij);
    eval_Pnv(t, f, ode, obs, Pnvij);
    eval_Pvv(t, f, ode, obs, Pvvij);    
    eval_Pnnn(t, f, ode, obs, Pnnnij);
    eval_Pnnv(t, f, ode, obs, Pnnvij);
    eval_Pnvv(t, f, ode, obs, Pnvvij);    
    eval_Pvvv(t, f, ode, obs, Pvvvij);    

    REAL hp22 = 0;
    hp22 += sqrt(PI/5.0)*(Rij[3*0+0] - Rij[3*1+1]);

    hp22 += ((2*sqrt(PI/5.0)*(nz*Pnij[3*0+0] - nx*Pnij[3*0+2] - nz*Pnij[3*1+1] + ny*Pnij[3*1+2] + vz*Pvij[3*0+0] - vx*Pvij[3*0+2] - vz*Pvij[3*1+1] + vy*Pvij[3*1+2]))/3.0);

   hp22 += ((sqrt(PI/5.0)*(3*nx*nx*Pnnij[3*0+0] + 7*ny*ny*Pnnij[3*0+0] + 11*nz*nz*Pnnij[3*0+0] - 12*nx*nz*Pnnij[3*0+2] - 7*nx*nx*Pnnij[3*1+1] - 3*ny*ny*Pnnij[3*1+1] - 11*nz*nz*Pnnij[3*1+1] + 12*ny*nz*Pnnij[3*1+2] + 4*nx*nx*Pnnij[3*2+2] - 4*ny*ny*Pnnij[3*2+2] + 3*nx*vx*Pnvij[3*0+0] + 7*ny*vy*Pnvij[3*0+0] + 11*nz*vz*Pnvij[3*0+0] - 6*nz*vx*Pnvij[3*0+2] - 6*nx*vz*Pnvij[3*0+2] - 7*nx*vx*Pnvij[3*1+1] - 3*ny*vy*Pnvij[3*1+1] - 11*nz*vz*Pnvij[3*1+1] + 6*nz*vy*Pnvij[3*1+2] + 6*ny*vz*Pnvij[3*1+2] + 4*nx*vx*Pnvij[3*2+2] - 4*ny*vy*Pnvij[3*2+2] + 3*vx*vx*Pvvij[3*0+0] + 7*vy*vy*Pvvij[3*0+0] + 11*vz*vz*Pvvij[3*0+0] - 12*vx*vz*Pvvij[3*0+2] - 7*vx*vx*Pvvij[3*1+1] - 3*vy*vy*Pvvij[3*1+1] - 11*vz*vz*Pvvij[3*1+1] + 12*vy*vz*Pvvij[3*1+2] + 4*vx*vx*Pvvij[3*2+2] - 4*vy*vy*Pvvij[3*2+2]))/21.0);

    hp22 += (-(sqrt(PI/5.0)*(-6*nx*nx*nz*Pnnnij[3*0+0] - 
     9*ny*ny*nz*Pnnnij[3*0+0] - 9*nz*nz*nz*Pnnnij[3*0+0] + 
     6*nx*nx*nx*Pnnnij[3*0+2] + 12*nx*ny*ny*Pnnnij[3*0+2] + 
     12*nx*nz*nz*Pnnnij[3*0+2] + 9*nx*nx*nz*Pnnnij[3*1+1] + 
     6*ny*ny*nz*Pnnnij[3*1+1] + 9*nz*nz*nz*Pnnnij[3*1+1] - 
     12*nx*nx*ny*Pnnnij[3*1+2] - 6*ny*ny*ny*Pnnnij[3*1+2] - 
     12*ny*nz*nz*Pnnnij[3*1+2] - 3*nx*nx*nz*Pnnnij[3*2+2] + 
     3*ny*ny*nz*Pnnnij[3*2+2] - 4*nx*nz*vx*Pnnvij[3*0+0] - 
     6*ny*nz*vy*Pnnvij[3*0+0] - 2*nx*nx*vz*Pnnvij[3*0+0] - 
     3*ny*ny*vz*Pnnvij[3*0+0] - 9*nz*nz*vz*Pnnvij[3*0+0] + 
     6*nx*nx*vx*Pnnvij[3*0+2] + 4*ny*ny*vx*Pnnvij[3*0+2] + 
     4*nz*nz*vx*Pnnvij[3*0+2] + 8*nx*ny*vy*Pnnvij[3*0+2] + 
     8*nx*nz*vz*Pnnvij[3*0+2] + 6*nx*nz*vx*Pnnvij[3*1+1] + 
     4*ny*nz*vy*Pnnvij[3*1+1] + 3*nx*nx*vz*Pnnvij[3*1+1] + 
     2*ny*ny*vz*Pnnvij[3*1+1] + 9*nz*nz*vz*Pnnvij[3*1+1] - 
     8*nx*ny*vx*Pnnvij[3*1+2] - 4*nx*nx*vy*Pnnvij[3*1+2] - 
     6*ny*ny*vy*Pnnvij[3*1+2] - 4*nz*nz*vy*Pnnvij[3*1+2] - 
     8*ny*nz*vz*Pnnvij[3*1+2] - 2*nx*nz*vx*Pnnvij[3*2+2] + 
     2*ny*nz*vy*Pnnvij[3*2+2] - nx*nx*vz*Pnnvij[3*2+2] + 
     ny*ny*vz*Pnnvij[3*2+2] - 2*nz*vx*vx*Pnvvij[3*0+0] - 
     3*nz*vy*vy*Pnvvij[3*0+0] - 4*nx*vx*vz*Pnvvij[3*0+0] - 
     6*ny*vy*vz*Pnvvij[3*0+0] - 9*nz*vz*vz*Pnvvij[3*0+0] + 
     6*nx*vx*vx*Pnvvij[3*0+2] + 8*ny*vx*vy*Pnvvij[3*0+2] + 
     4*nx*vy*vy*Pnvvij[3*0+2] + 8*nz*vx*vz*Pnvvij[3*0+2] + 
     4*nx*vz*vz*Pnvvij[3*0+2] + 3*nz*vx*vx*Pnvvij[3*1+1] + 
     2*nz*vy*vy*Pnvvij[3*1+1] + 6*nx*vx*vz*Pnvvij[3*1+1] + 
     4*ny*vy*vz*Pnvvij[3*1+1] + 9*nz*vz*vz*Pnvvij[3*1+1] - 
     4*ny*vx*vx*Pnvvij[3*1+2] - 8*nx*vx*vy*Pnvvij[3*1+2] - 
     6*ny*vy*vy*Pnvvij[3*1+2] - 8*nz*vy*vz*Pnvvij[3*1+2] - 
     4*ny*vz*vz*Pnvvij[3*1+2] - nz*vx*vx*Pnvvij[3*2+2] + 
     nz*vy*vy*Pnvvij[3*2+2] - 2*nx*vx*vz*Pnvvij[3*2+2] + 
     2*ny*vy*vz*Pnvvij[3*2+2] - 6*vx*vx*vz*Pvvvij[3*0+0] - 
     9*vy*vy*vz*Pvvvij[3*0+0] - 9*vz*vz*vz*Pvvvij[3*0+0] + 
     6*vx*vx*vx*Pvvvij[3*0+2] + 12*vx*vy*vy*Pvvvij[3*0+2] + 
     12*vx*vz*vz*Pvvvij[3*0+2] + 9*vx*vx*vz*Pvvvij[3*1+1] + 
     6*vy*vy*vz*Pvvvij[3*1+1] + 9*vz*vz*vz*Pvvvij[3*1+1] - 
     12*vx*vx*vy*Pvvvij[3*1+2] - 6*vy*vy*vy*Pvvvij[3*1+2] - 
     12*vy*vz*vz*Pvvvij[3*1+2] - 3*vx*vx*vz*Pvvvij[3*2+2] + 
     3*vy*vy*vz*Pvvvij[3*2+2]))/21.0);

    if((ode.hterms & H_PQSO) != 0) {
     hp22 += 4*ode.mu/rsq*sqrt(PI/5.0)*(nx * Dy + ny * Dx);
    }

    if((ode.hterms & H_P15QSO) != 0) {
     hp22 += 2*ode.mu*((2*sqrt(PI/5.0)*(-6*dm*Dz*nx*ny*rdot + 3*dm*Dy*nx*nz*rdot + 
      3*dm*Dx*ny*nz*rdot + 3*ode.m*ny*nz*rdot*Sx + 3*ode.m*nx*nz*rdot*Sy - 
      6*ode.m*nx*ny*rdot*Sz + 2*dm*Dz*ny*vx + 2*dm*Dy*nz*vx + 2*ode.m*nz*Sy*vx + 
      2*ode.m*ny*Sz*vx + 2*dm*Dz*nx*vy + 2*dm*Dx*nz*vy + 2*ode.m*nz*Sx*vy + 
      2*ode.m*nx*Sz*vy - 4*dm*Dy*nx*vz - 4*dm*Dx*ny*vz - 4*ode.m*ny*Sx*vz - 
      4*ode.m*nx*Sy*vz))/(3.0*ode.m*r*r));
     hp22 += 2*ode.mu*((4*sqrt(PI/5.0)*(dm*Dz*ny*vx - 2*dm*Dy*nz*vx - 2*ode.m*nz*Sy*vx + 
      ode.m*ny*Sz*vx + dm*Dz*nx*vy - 2*dm*Dx*nz*vy - 2*ode.m*nz*Sx*vy + ode.m*nx*Sz*vy + 
      dm*Dy*nx*vz + dm*Dx*ny*vz + ode.m*ny*Sx*vz + ode.m*nx*Sy*vz))/(3.0*ode.m*r*r));
    }
   return hp22;
}

static FUNCDEF(hp2m2) {
    REAL Rij[9];
    REAL Pnij[9];
    REAL Pvij[9];
    REAL Pnnij[9];
    REAL Pnvij[9];
    REAL Pvvij[9];
    REAL Pnnnij[9];
    REAL Pnnvij[9];
    REAL Pnvvij[9];
    REAL Pvvvij[9];
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);
    REAL nx = rx/r;  
    REAL ny = ry/r;  
    REAL nz = rz/r;
   
    REAL vx = eval_vx(t, f, ode, obs);
    REAL vy = eval_vy(t, f, ode, obs);    
    REAL vz = eval_vz(t, f, ode, obs);  

    REAL rdot = (rx*vx + ry*vy + rz*vz)/r;
    REAL dm = ode.m1 - ode.m2;

    REAL cSs1 = ode.m1*ode.m1;
    REAL S1x = cSs1*f[I_s1x];
    REAL S1y = cSs1*f[I_s1y];
    REAL S1z = cSs1*f[I_s1z];
    REAL cSs2 = ode.m2*ode.m2;
    REAL S2x = cSs2*f[I_s2x];
    REAL S2y = cSs2*f[I_s2y];
    REAL S2z = cSs2*f[I_s2z];
    REAL Sx = S1x + S2x;
    REAL Sy = S1y + S2y;
    REAL Sz = S1z + S2z;
    REAL Dx = ode.m*(S2x/ode.m2 - S1x/ode.m1);
    REAL Dy = ode.m*(S2y/ode.m2 - S1y/ode.m1);
    REAL Dz = ode.m*(S2z/ode.m2 - S1z/ode.m1);

    eval_R(t, f, ode, obs, Rij);
    eval_Pn(t, f, ode, obs, Pnij);
    eval_Pv(t, f, ode, obs, Pvij);    
    eval_Pnn(t, f, ode, obs, Pnnij);
    eval_Pnv(t, f, ode, obs, Pnvij);
    eval_Pvv(t, f, ode, obs, Pvvij);    
    eval_Pnnn(t, f, ode, obs, Pnnnij);
    eval_Pnnv(t, f, ode, obs, Pnnvij);
    eval_Pnvv(t, f, ode, obs, Pnvvij);    
    eval_Pvvv(t, f, ode, obs, Pvvvij);    

    REAL hp2m2 = 0;
    hp2m2 += sqrt(PI/5.0)*(Rij[3*0+0] - Rij[3*1+1]);

    hp2m2 += ((-2*sqrt(PI/5.0)*(nz*Pnij[3*0+0] - nx*Pnij[3*0+2] - nz*Pnij[3*1+1] + ny*Pnij[3*1+2] + vz*Pvij[3*0+0] - vx*Pvij[3*0+2] - vz*Pvij[3*1+1] + vy*Pvij[3*1+2]))/3.0);

    hp2m2 += ((sqrt(PI/5.0)*(3*nx*nx*Pnnij[3*0+0] + 7*ny*ny*Pnnij[3*0+0] + 11*nz*nz*Pnnij[3*0+0] - 12*nx*nz*Pnnij[3*0+2] - 7*nx*nx*Pnnij[3*1+1] - 3*ny*ny*Pnnij[3*1+1] - 11*nz*nz*Pnnij[3*1+1] + 12*ny*nz*Pnnij[3*1+2] + 4*nx*nx*Pnnij[3*2+2] - 4*ny*ny*Pnnij[3*2+2] + 3*nx*vx*Pnvij[3*0+0] + 7*ny*vy*Pnvij[3*0+0] + 11*nz*vz*Pnvij[3*0+0] - 6*nz*vx*Pnvij[3*0+2] - 6*nx*vz*Pnvij[3*0+2] - 7*nx*vx*Pnvij[3*1+1] - 3*ny*vy*Pnvij[3*1+1] - 11*nz*vz*Pnvij[3*1+1] + 6*nz*vy*Pnvij[3*1+2] + 6*ny*vz*Pnvij[3*1+2] + 4*nx*vx*Pnvij[3*2+2] - 4*ny*vy*Pnvij[3*2+2] + 3*vx*vx*Pvvij[3*0+0] + 7*vy*vy*Pvvij[3*0+0] + 11*vz*vz*Pvvij[3*0+0] - 12*vx*vz*Pvvij[3*0+2] - 7*vx*vx*Pvvij[3*1+1] - 3*vy*vy*Pvvij[3*1+1] - 11*vz*vz*Pvvij[3*1+1] + 12*vy*vz*Pvvij[3*1+2] + 4*vx*vx*Pvvij[3*2+2] - 4*vy*vy*Pvvij[3*2+2]))/21.0);

    hp2m2 += (-(sqrt(PI/5.0)*(6*nx*nx*nz*Pnnnij[3*0+0] + 
     9*ny*ny*nz*Pnnnij[3*0+0] + 9*nz*nz*nz*Pnnnij[3*0+0] - 
     6*nx*nx*nx*Pnnnij[3*0+2] - 12*nx*ny*ny*Pnnnij[3*0+2] - 
     12*nx*nz*nz*Pnnnij[3*0+2] - 9*nx*nx*nz*Pnnnij[3*1+1] - 
     6*ny*ny*nz*Pnnnij[3*1+1] - 9*nz*nz*nz*Pnnnij[3*1+1] + 
     12*nx*nx*ny*Pnnnij[3*1+2] + 6*ny*ny*ny*Pnnnij[3*1+2] + 
     12*ny*nz*nz*Pnnnij[3*1+2] + 3*nx*nx*nz*Pnnnij[3*2+2] - 
     3*ny*ny*nz*Pnnnij[3*2+2] + 4*nx*nz*vx*Pnnvij[3*0+0] + 
     6*ny*nz*vy*Pnnvij[3*0+0] + 2*nx*nx*vz*Pnnvij[3*0+0] + 
     3*ny*ny*vz*Pnnvij[3*0+0] + 9*nz*nz*vz*Pnnvij[3*0+0] - 
     6*nx*nx*vx*Pnnvij[3*0+2] - 4*ny*ny*vx*Pnnvij[3*0+2] - 
     4*nz*nz*vx*Pnnvij[3*0+2] - 8*nx*ny*vy*Pnnvij[3*0+2] - 
     8*nx*nz*vz*Pnnvij[3*0+2] - 6*nx*nz*vx*Pnnvij[3*1+1] - 
     4*ny*nz*vy*Pnnvij[3*1+1] - 3*nx*nx*vz*Pnnvij[3*1+1] - 
     2*ny*ny*vz*Pnnvij[3*1+1] - 9*nz*nz*vz*Pnnvij[3*1+1] + 
     8*nx*ny*vx*Pnnvij[3*1+2] + 4*nx*nx*vy*Pnnvij[3*1+2] + 
     6*ny*ny*vy*Pnnvij[3*1+2] + 4*nz*nz*vy*Pnnvij[3*1+2] + 
     8*ny*nz*vz*Pnnvij[3*1+2] + 2*nx*nz*vx*Pnnvij[3*2+2] - 
     2*ny*nz*vy*Pnnvij[3*2+2] + nx*nx*vz*Pnnvij[3*2+2] - 
     ny*ny*vz*Pnnvij[3*2+2] + 2*nz*vx*vx*Pnvvij[3*0+0] + 
     3*nz*vy*vy*Pnvvij[3*0+0] + 4*nx*vx*vz*Pnvvij[3*0+0] + 
     6*ny*vy*vz*Pnvvij[3*0+0] + 9*nz*vz*vz*Pnvvij[3*0+0] - 
     6*nx*vx*vx*Pnvvij[3*0+2] - 8*ny*vx*vy*Pnvvij[3*0+2] - 
     4*nx*vy*vy*Pnvvij[3*0+2] - 8*nz*vx*vz*Pnvvij[3*0+2] - 
     4*nx*vz*vz*Pnvvij[3*0+2] - 3*nz*vx*vx*Pnvvij[3*1+1] - 
     2*nz*vy*vy*Pnvvij[3*1+1] - 6*nx*vx*vz*Pnvvij[3*1+1] - 
     4*ny*vy*vz*Pnvvij[3*1+1] - 9*nz*vz*vz*Pnvvij[3*1+1] + 
     4*ny*vx*vx*Pnvvij[3*1+2] + 8*nx*vx*vy*Pnvvij[3*1+2] + 
     6*ny*vy*vy*Pnvvij[3*1+2] + 8*nz*vy*vz*Pnvvij[3*1+2] + 
     4*ny*vz*vz*Pnvvij[3*1+2] + nz*vx*vx*Pnvvij[3*2+2] - 
     nz*vy*vy*Pnvvij[3*2+2] + 2*nx*vx*vz*Pnvvij[3*2+2] - 
     2*ny*vy*vz*Pnvvij[3*2+2] + 6*vx*vx*vz*Pvvvij[3*0+0] + 
     9*vy*vy*vz*Pvvvij[3*0+0] + 9*vz*vz*vz*Pvvvij[3*0+0] - 
     6*vx*vx*vx*Pvvvij[3*0+2] - 12*vx*vy*vy*Pvvvij[3*0+2] - 
     12*vx*vz*vz*Pvvvij[3*0+2] - 9*vx*vx*vz*Pvvvij[3*1+1] - 
     6*vy*vy*vz*Pvvvij[3*1+1] - 9*vz*vz*vz*Pvvvij[3*1+1] + 
     12*vx*vx*vy*Pvvvij[3*1+2] + 6*vy*vy*vy*Pvvvij[3*1+2] + 
     12*vy*vz*vz*Pvvvij[3*1+2] + 3*vx*vx*vz*Pvvvij[3*2+2] - 
     3*vy*vy*vz*Pvvvij[3*2+2]))/21.0);

    if((ode.hterms & H_PQSO) != 0) {
     hp2m2 += -4*ode.mu/rsq*sqrt(PI/5.0)*(nx * Dy + ny * Dx);
    }

    if((ode.hterms & H_P15QSO) != 0) {
     hp2m2 += 2*ode.mu*((2*sqrt(PI/5.0)*(-6*dm*Dz*nx*ny*rdot + 3*dm*Dy*nx*nz*rdot + 
    3*dm*Dx*ny*nz*rdot + 3*ode.m*ny*nz*rdot*Sx + 3*ode.m*nx*nz*rdot*Sy - 
    6*ode.m*nx*ny*rdot*Sz + 2*dm*Dz*ny*vx + 2*dm*Dy*nz*vx + 2*ode.m*nz*Sy*vx + 
    2*ode.m*ny*Sz*vx + 2*dm*Dz*nx*vy + 2*dm*Dx*nz*vy + 2*ode.m*nz*Sx*vy + 
    2*ode.m*nx*Sz*vy - 4*dm*Dy*nx*vz - 4*dm*Dx*ny*vz - 4*ode.m*ny*Sx*vz - 
    4*ode.m*nx*Sy*vz))/(3.0*ode.m*r*r));
     hp2m2 += 2*ode.mu*((4*sqrt(PI/5.0)*(dm*Dz*ny*vx - 2*dm*Dy*nz*vx - 2*ode.m*nz*Sy*vx + 
    ode.m*ny*Sz*vx + dm*Dz*nx*vy - 2*dm*Dx*nz*vy - 2*ode.m*nz*Sx*vy + ode.m*nx*Sz*vy + 
    dm*Dy*nx*vz + dm*Dx*ny*vz + ode.m*ny*Sx*vz + ode.m*nx*Sy*vz))/(3.0*ode.m*r*r));
    }
   return hp2m2;
}

static FUNCDEF(hx22) {
    REAL Rij[9];
    REAL Pnij[9];
    REAL Pvij[9];
    REAL Pnnij[9];
    REAL Pnvij[9];
    REAL Pvvij[9];
    REAL Pnnnij[9];
    REAL Pnnvij[9];
    REAL Pnvvij[9];
    REAL Pvvvij[9];
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);
    REAL nx = rx/r;  
    REAL ny = ry/r;  
    REAL nz = rz/r;
   
    REAL vx = eval_vx(t, f, ode, obs);
    REAL vy = eval_vy(t, f, ode, obs);    
    REAL vz = eval_vz(t, f, ode, obs);

    REAL rdot = (rx*vx + ry*vy + rz*vz)/r;
    REAL dm = ode.m1 - ode.m2;

    REAL cSs1 = ode.m1*ode.m1;
    REAL S1x = cSs1*f[I_s1x];
    REAL S1y = cSs1*f[I_s1y];
    REAL S1z = cSs1*f[I_s1z];
    REAL cSs2 = ode.m2*ode.m2;
    REAL S2x = cSs2*f[I_s2x];
    REAL S2y = cSs2*f[I_s2y];
    REAL S2z = cSs2*f[I_s2z];
    REAL Sx = S1x + S2x;
    REAL Sy = S1y + S2y;
    REAL Sz = S1z + S2z;
    REAL Dx = ode.m*(S2x/ode.m2 - S1x/ode.m1);
    REAL Dy = ode.m*(S2y/ode.m2 - S1y/ode.m1);
    REAL Dz = ode.m*(S2z/ode.m2 - S1z/ode.m1);

    eval_R(t, f, ode, obs, Rij);
    eval_Pn(t, f, ode, obs, Pnij);    
    eval_Pv(t, f, ode, obs, Pvij);
    eval_Pnn(t, f, ode, obs, Pnnij);    
    eval_Pnv(t, f, ode, obs, Pnvij);    
    eval_Pvv(t, f, ode, obs, Pvvij);
    eval_Pnnn(t, f, ode, obs, Pnnnij);
    eval_Pnnv(t, f, ode, obs, Pnnvij);
    eval_Pnvv(t, f, ode, obs, Pnvvij);    
    eval_Pvvv(t, f, ode, obs, Pvvvij);    

    REAL hx22 = 0;  
    hx22 += 2*sqrt(PI/5.0)*Rij[3*0+1];

    hx22 += ((2*sqrt(PI/5.0)*(2*nz*Pnij[3*0+1] - ny*Pnij[3*0+2] - nx*Pnij[3*1+2] + 2*vz*Pvij[3*0+1] - vy*Pvij[3*0+2] - vx*Pvij[3*1+2]))/3.0);

    hx22 += ((2*sqrt(PI/5.0)*(-2*nx*ny*Pnnij[3*0+0] + 5*nx*nx*Pnnij[3*0+1] + 5*ny*ny*Pnnij[3*0+1] + 11*nz*nz*Pnnij[3*0+1] - 6*ny*nz*Pnnij[3*0+2] - 2*nx*ny*Pnnij[3*1+1] - 6*nx*nz*Pnnij[3*1+2] + 4*nx*ny*Pnnij[3*2+2] - ny*vx*Pnvij[3*0+0] - nx*vy*Pnvij[3*0+0] + 5*nx*vx*Pnvij[3*0+1] + 5*ny*vy*Pnvij[3*0+1] + 11*nz*vz*Pnvij[3*0+1] - 3*nz*vy*Pnvij[3*0+2] - 3*ny*vz*Pnvij[3*0+2] - ny*vx*Pnvij[3*1+1] - nx*vy*Pnvij[3*1+1] - 3*nz*vx*Pnvij[3*1+2] - 3*nx*vz*Pnvij[3*1+2] + 2*ny*vx*Pnvij[3*2+2] + 2*nx*vy*Pnvij[3*2+2] - 2*vx*vy*Pvvij[3*0+0] + 5*vx*vx*Pvvij[3*0+1] + 5*vy*vy*Pvvij[3*0+1] + 11*vz*vz*Pvvij[3*0+1] - 6*vy*vz*Pvvij[3*0+2] - 2*vx*vy*Pvvij[3*1+1] - 6*vx*vz*Pvvij[3*1+2] + 4*vx*vy*Pvvij[3*2+2]))/21.0);
    hx22 += (-(sqrt(PI/5.0)*(3*nx*ny*nz*Pnnnij[3*0+0] - 
     15*nx*nx*nz*Pnnnij[3*0+1] - 15*ny*ny*nz*Pnnnij[3*0+1] - 
     18*nz*nz*nz*Pnnnij[3*0+1] + 3*nx*nx*ny*Pnnnij[3*0+2] + 
     9*ny*ny*ny*Pnnnij[3*0+2] + 12*ny*nz*nz*Pnnnij[3*0+2] + 
     3*nx*ny*nz*Pnnnij[3*1+1] + 9*nx*nx*nx*Pnnnij[3*1+2] + 
     3*nx*ny*ny*Pnnnij[3*1+2] + 12*nx*nz*nz*Pnnnij[3*1+2] - 
     6*nx*ny*nz*Pnnnij[3*2+2] + ny*nz*vx*Pnnvij[3*0+0] + 
     nx*nz*vy*Pnnvij[3*0+0] + nx*ny*vz*Pnnvij[3*0+0] - 
     10*nx*nz*vx*Pnnvij[3*0+1] - 10*ny*nz*vy*Pnnvij[3*0+1] - 
     5*nx*nx*vz*Pnnvij[3*0+1] - 5*ny*ny*vz*Pnnvij[3*0+1] - 
     18*nz*nz*vz*Pnnvij[3*0+1] + 2*nx*ny*vx*Pnnvij[3*0+2] + 
     nx*nx*vy*Pnnvij[3*0+2] + 9*ny*ny*vy*Pnnvij[3*0+2] + 
     4*nz*nz*vy*Pnnvij[3*0+2] + 8*ny*nz*vz*Pnnvij[3*0+2] + 
     ny*nz*vx*Pnnvij[3*1+1] + nx*nz*vy*Pnnvij[3*1+1] + 
     nx*ny*vz*Pnnvij[3*1+1] + 9*nx*nx*vx*Pnnvij[3*1+2] + 
     ny*ny*vx*Pnnvij[3*1+2] + 4*nz*nz*vx*Pnnvij[3*1+2] + 
     2*nx*ny*vy*Pnnvij[3*1+2] + 8*nx*nz*vz*Pnnvij[3*1+2] - 
     2*ny*nz*vx*Pnnvij[3*2+2] - 2*nx*nz*vy*Pnnvij[3*2+2] - 
     2*nx*ny*vz*Pnnvij[3*2+2] + nz*vx*vy*Pnvvij[3*0+0] + 
     ny*vx*vz*Pnvvij[3*0+0] + nx*vy*vz*Pnvvij[3*0+0] - 
     5*nz*vx*vx*Pnvvij[3*0+1] - 5*nz*vy*vy*Pnvvij[3*0+1] - 
     10*nx*vx*vz*Pnvvij[3*0+1] - 10*ny*vy*vz*Pnvvij[3*0+1] - 
     18*nz*vz*vz*Pnvvij[3*0+1] + ny*vx*vx*Pnvvij[3*0+2] + 
     2*nx*vx*vy*Pnvvij[3*0+2] + 9*ny*vy*vy*Pnvvij[3*0+2] + 
     8*nz*vy*vz*Pnvvij[3*0+2] + 4*ny*vz*vz*Pnvvij[3*0+2] + 
     nz*vx*vy*Pnvvij[3*1+1] + ny*vx*vz*Pnvvij[3*1+1] + 
     nx*vy*vz*Pnvvij[3*1+1] + 9*nx*vx*vx*Pnvvij[3*1+2] + 
     2*ny*vx*vy*Pnvvij[3*1+2] + nx*vy*vy*Pnvvij[3*1+2] + 
     8*nz*vx*vz*Pnvvij[3*1+2] + 4*nx*vz*vz*Pnvvij[3*1+2] - 
     2*nz*vx*vy*Pnvvij[3*2+2] - 2*ny*vx*vz*Pnvvij[3*2+2] - 
     2*nx*vy*vz*Pnvvij[3*2+2] + 3*vx*vy*vz*Pvvvij[3*0+0] - 
     15*vx*vx*vz*Pvvvij[3*0+1] - 15*vy*vy*vz*Pvvvij[3*0+1] - 
     18*vz*vz*vz*Pvvvij[3*0+1] + 3*vx*vx*vy*Pvvvij[3*0+2] + 
     9*vy*vy*vy*Pvvvij[3*0+2] + 12*vy*vz*vz*Pvvvij[3*0+2] + 
     3*vx*vy*vz*Pvvvij[3*1+1] + 9*vx*vx*vx*Pvvvij[3*1+2] + 
     3*vx*vy*vy*Pvvvij[3*1+2] + 12*vx*vz*vz*Pvvvij[3*1+2] - 
     6*vx*vy*vz*Pvvvij[3*2+2]))/21.0);

    if((ode.hterms & H_PQSO) != 0) {
     hx22 += -4*ode.mu/rsq*sqrt(PI/5.0)*(nx * Dx - ny * Dy);
    }

    if((ode.hterms & H_P15QSO) != 0) {
     hx22 += 2*ode.mu*((2*sqrt(PI/5.0)*(3*dm*Dz*nx*nx*rdot - 3*dm*Dz*ny*ny*rdot - 
    3*dm*Dx*nx*nz*rdot + 3*dm*Dy*ny*nz*rdot - 3*ode.m*nx*nz*rdot*Sx + 
    3*ode.m*ny*nz*rdot*Sy + 3*ode.m*nx*nx*rdot*Sz - 3*ode.m*ny*ny*rdot*Sz - 2*dm*Dz*nx*vx - 
    2*dm*Dx*nz*vx - 2*ode.m*nz*Sx*vx - 2*ode.m*nx*Sz*vx + 2*dm*Dz*ny*vy + 
    2*dm*Dy*nz*vy + 2*ode.m*nz*Sy*vy + 2*ode.m*ny*Sz*vy + 4*dm*Dx*nx*vz - 
    4*dm*Dy*ny*vz + 4*ode.m*nx*Sx*vz - 4*ode.m*ny*Sy*vz))/(3.0*ode.m*r*r));
     hx22 += 2*ode.mu*((-4*sqrt(PI/5.0)*(dm*Dz*nx*vx - 2*dm*Dx*nz*vx - 2*ode.m*nz*Sx*vx + 
    ode.m*nx*Sz*vx - dm*Dz*ny*vy + 2*dm*Dy*nz*vy + 2*ode.m*nz*Sy*vy - ode.m*ny*Sz*vy + 
    dm*Dx*nx*vz - dm*Dy*ny*vz + ode.m*nx*Sx*vz - ode.m*ny*Sy*vz))/(3.0*ode.m*r*r));
    }
    return hx22;
}

static FUNCDEF(hx2m2) {
    REAL Rij[9];
    REAL Pnij[9];
    REAL Pvij[9];
    REAL Pnnij[9];
    REAL Pnvij[9];
    REAL Pvvij[9];
    REAL Pnnnij[9];
    REAL Pnnvij[9];
    REAL Pnvvij[9];
    REAL Pvvvij[9];
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);
    REAL nx = rx/r;  
    REAL ny = ry/r;  
    REAL nz = rz/r;
   
    REAL vx = eval_vx(t, f, ode, obs);
    REAL vy = eval_vy(t, f, ode, obs);    
    REAL vz = eval_vz(t, f, ode, obs);

    REAL rdot = (rx*vx + ry*vy + rz*vz)/r;
    REAL dm = ode.m1 - ode.m2;

    REAL cSs1 = ode.m1*ode.m1;
    REAL S1x = cSs1*f[I_s1x];
    REAL S1y = cSs1*f[I_s1y];
    REAL S1z = cSs1*f[I_s1z];
    REAL cSs2 = ode.m2*ode.m2;
    REAL S2x = cSs2*f[I_s2x];
    REAL S2y = cSs2*f[I_s2y];
    REAL S2z = cSs2*f[I_s2z];
    REAL Sx = S1x + S2x;
    REAL Sy = S1y + S2y;
    REAL Sz = S1z + S2z;
    REAL Dx = ode.m*(S2x/ode.m2 - S1x/ode.m1);
    REAL Dy = ode.m*(S2y/ode.m2 - S1y/ode.m1);
    REAL Dz = ode.m*(S2z/ode.m2 - S1z/ode.m1);

    eval_R(t, f, ode, obs, Rij);
    eval_Pn(t, f, ode, obs, Pnij);    
    eval_Pv(t, f, ode, obs, Pvij);
    eval_Pnn(t, f, ode, obs, Pnnij);    
    eval_Pnv(t, f, ode, obs, Pnvij);    
    eval_Pvv(t, f, ode, obs, Pvvij);
    eval_Pnnn(t, f, ode, obs, Pnnnij);
    eval_Pnnv(t, f, ode, obs, Pnnvij);
    eval_Pnvv(t, f, ode, obs, Pnvvij);    
    eval_Pvvv(t, f, ode, obs, Pvvvij);    

    REAL hx2m2 = 0;  
    hx2m2 += -2*sqrt(PI/5.0)*Rij[3*0+1];

    hx2m2 += ((2*sqrt(PI/5.0)*(2*nz*Pnij[3*0+1] - ny*Pnij[3*0+2] - nx*Pnij[3*1+2] + 2*vz*Pvij[3*0+1] - vy*Pvij[3*0+2] - vx*Pvij[3*1+2]))/3.0);

    hx2m2 += ((-2*sqrt(PI/5.0)*(-2*nx*ny*Pnnij[3*0+0] + 5*nx*nx*Pnnij[3*0+1] + 5*ny*ny*Pnnij[3*0+1] + 11*nz*nz*Pnnij[3*0+1] - 6*ny*nz*Pnnij[3*0+2] - 2*nx*ny*Pnnij[3*1+1] - 6*nx*nz*Pnnij[3*1+2] + 4*nx*ny*Pnnij[3*2+2] - ny*vx*Pnvij[3*0+0] - nx*vy*Pnvij[3*0+0] + 5*nx*vx*Pnvij[3*0+1] + 5*ny*vy*Pnvij[3*0+1] + 11*nz*vz*Pnvij[3*0+1] - 3*nz*vy*Pnvij[3*0+2] - 3*ny*vz*Pnvij[3*0+2] - ny*vx*Pnvij[3*1+1] - nx*vy*Pnvij[3*1+1] - 3*nz*vx*Pnvij[3*1+2] - 3*nx*vz*Pnvij[3*1+2] + 2*ny*vx*Pnvij[3*2+2] + 2*nx*vy*Pnvij[3*2+2] - 2*vx*vy*Pvvij[3*0+0] + 5*vx*vx*Pvvij[3*0+1] + 5*vy*vy*Pvvij[3*0+1] + 11*vz*vz*Pvvij[3*0+1] - 6*vy*vz*Pvvij[3*0+2] - 2*vx*vy*Pvvij[3*1+1] - 6*vx*vz*Pvvij[3*1+2] + 4*vx*vy*Pvvij[3*2+2]))/21.0);

    hx2m2 += (-(sqrt(PI/5.0)*(3*nx*ny*nz*Pnnnij[3*0+0] - 
     15*nx*nx*nz*Pnnnij[3*0+1] - 15*ny*ny*nz*Pnnnij[3*0+1] - 
     18*nz*nz*nz*Pnnnij[3*0+1] + 3*nx*nx*ny*Pnnnij[3*0+2] + 
     9*ny*ny*ny*Pnnnij[3*0+2] + 12*ny*nz*nz*Pnnnij[3*0+2] + 
     3*nx*ny*nz*Pnnnij[3*1+1] + 9*nx*nx*nx*Pnnnij[3*1+2] + 
     3*nx*ny*ny*Pnnnij[3*1+2] + 12*nx*nz*nz*Pnnnij[3*1+2] - 
     6*nx*ny*nz*Pnnnij[3*2+2] + ny*nz*vx*Pnnvij[3*0+0] + 
     nx*nz*vy*Pnnvij[3*0+0] + nx*ny*vz*Pnnvij[3*0+0] - 
     10*nx*nz*vx*Pnnvij[3*0+1] - 10*ny*nz*vy*Pnnvij[3*0+1] - 
     5*nx*nx*vz*Pnnvij[3*0+1] - 5*ny*ny*vz*Pnnvij[3*0+1] - 
     18*nz*nz*vz*Pnnvij[3*0+1] + 2*nx*ny*vx*Pnnvij[3*0+2] + 
     nx*nx*vy*Pnnvij[3*0+2] + 9*ny*ny*vy*Pnnvij[3*0+2] + 
     4*nz*nz*vy*Pnnvij[3*0+2] + 8*ny*nz*vz*Pnnvij[3*0+2] + 
     ny*nz*vx*Pnnvij[3*1+1] + nx*nz*vy*Pnnvij[3*1+1] + 
     nx*ny*vz*Pnnvij[3*1+1] + 9*nx*nx*vx*Pnnvij[3*1+2] + 
     ny*ny*vx*Pnnvij[3*1+2] + 4*nz*nz*vx*Pnnvij[3*1+2] + 
     2*nx*ny*vy*Pnnvij[3*1+2] + 8*nx*nz*vz*Pnnvij[3*1+2] - 
     2*ny*nz*vx*Pnnvij[3*2+2] - 2*nx*nz*vy*Pnnvij[3*2+2] - 
     2*nx*ny*vz*Pnnvij[3*2+2] + nz*vx*vy*Pnvvij[3*0+0] + 
     ny*vx*vz*Pnvvij[3*0+0] + nx*vy*vz*Pnvvij[3*0+0] - 
     5*nz*vx*vx*Pnvvij[3*0+1] - 5*nz*vy*vy*Pnvvij[3*0+1] - 
     10*nx*vx*vz*Pnvvij[3*0+1] - 10*ny*vy*vz*Pnvvij[3*0+1] - 
     18*nz*vz*vz*Pnvvij[3*0+1] + ny*vx*vx*Pnvvij[3*0+2] + 
     2*nx*vx*vy*Pnvvij[3*0+2] + 9*ny*vy*vy*Pnvvij[3*0+2] + 
     8*nz*vy*vz*Pnvvij[3*0+2] + 4*ny*vz*vz*Pnvvij[3*0+2] + 
     nz*vx*vy*Pnvvij[3*1+1] + ny*vx*vz*Pnvvij[3*1+1] + 
     nx*vy*vz*Pnvvij[3*1+1] + 9*nx*vx*vx*Pnvvij[3*1+2] + 
     2*ny*vx*vy*Pnvvij[3*1+2] + nx*vy*vy*Pnvvij[3*1+2] + 
     8*nz*vx*vz*Pnvvij[3*1+2] + 4*nx*vz*vz*Pnvvij[3*1+2] - 
     2*nz*vx*vy*Pnvvij[3*2+2] - 2*ny*vx*vz*Pnvvij[3*2+2] - 
     2*nx*vy*vz*Pnvvij[3*2+2] + 3*vx*vy*vz*Pvvvij[3*0+0] - 
     15*vx*vx*vz*Pvvvij[3*0+1] - 15*vy*vy*vz*Pvvvij[3*0+1] - 
     18*vz*vz*vz*Pvvvij[3*0+1] + 3*vx*vx*vy*Pvvvij[3*0+2] + 
     9*vy*vy*vy*Pvvvij[3*0+2] + 12*vy*vz*vz*Pvvvij[3*0+2] + 
     3*vx*vy*vz*Pvvvij[3*1+1] + 9*vx*vx*vx*Pvvvij[3*1+2] + 
     3*vx*vy*vy*Pvvvij[3*1+2] + 12*vx*vz*vz*Pvvvij[3*1+2] - 
     6*vx*vy*vz*Pvvvij[3*2+2]))/21.0);

    if((ode.hterms & H_PQSO) != 0) {
     hx2m2 += -4*ode.mu/rsq*sqrt(PI/5.0)*(nx * Dx - ny * Dy);
    }

    if((ode.hterms & H_P15QSO) != 0) {
     hx2m2 += 2*ode.mu*((-2*sqrt(PI/5.0)*(3*dm*Dz*nx*nx*rdot - 3*dm*Dz*ny*ny*rdot - 
    3*dm*Dx*nx*nz*rdot + 3*dm*Dy*ny*nz*rdot - 3*ode.m*nx*nz*rdot*Sx + 
    3*ode.m*ny*nz*rdot*Sy + 3*ode.m*nx*nx*rdot*Sz - 3*ode.m*ny*ny*rdot*Sz - 2*dm*Dz*nx*vx - 
    2*dm*Dx*nz*vx - 2*ode.m*nz*Sx*vx - 2*ode.m*nx*Sz*vx + 2*dm*Dz*ny*vy + 
    2*dm*Dy*nz*vy + 2*ode.m*nz*Sy*vy + 2*ode.m*ny*Sz*vy + 4*dm*Dx*nx*vz - 
    4*dm*Dy*ny*vz + 4*ode.m*nx*Sx*vz - 4*ode.m*ny*Sy*vz))/(3.0*ode.m*r*r));
     hx2m2 += 2*ode.mu*((4*sqrt(PI/5.0)*(dm*Dz*nx*vx - 2*dm*Dx*nz*vx - 2*ode.m*nz*Sx*vx + 
    ode.m*nx*Sz*vx - dm*Dz*ny*vy + 2*dm*Dy*nz*vy + 2*ode.m*nz*Sy*vy - ode.m*ny*Sz*vy + 
    dm*Dx*nx*vz - dm*Dy*ny*vz + ode.m*nx*Sx*vz - ode.m*ny*Sy*vz))/(3.0*ode.m*r*r));
    }
    return hx2m2;
}

static FUNCDEF(hp21) {
    REAL Rij[9];
    REAL Pnij[9];
    REAL Pvij[9];
    REAL Pnnij[9];
    REAL Pnvij[9];
    REAL Pvvij[9];
    REAL Pnnnij[9];
    REAL Pnnvij[9];
    REAL Pnvvij[9];
    REAL Pvvvij[9];
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);
    REAL nx = rx/r;  
    REAL ny = ry/r;  
    REAL nz = rz/r;

    REAL vx = eval_vx(t, f, ode, obs);
    REAL vy = eval_vy(t, f, ode, obs);    
    REAL vz = eval_vz(t, f, ode, obs);  

    REAL rdot = (rx*vx + ry*vy + rz*vz)/r;
    REAL dm = ode.m1 - ode.m2;
   
    REAL cSs1 = ode.m1*ode.m1;
    REAL S1x = cSs1*f[I_s1x];
    REAL S1y = cSs1*f[I_s1y];
    REAL S1z = cSs1*f[I_s1z];
    REAL cSs2 = ode.m2*ode.m2;
    REAL S2x = cSs2*f[I_s2x];
    REAL S2y = cSs2*f[I_s2y];
    REAL S2z = cSs2*f[I_s2z];
    REAL Sx = S1x + S2x;
    REAL Sy = S1y + S2y;
    REAL Sz = S1z + S2z;
    REAL Dx = ode.m*(S2x/ode.m2 - S1x/ode.m1);
    REAL Dy = ode.m*(S2y/ode.m2 - S1y/ode.m1);
    REAL Dz = ode.m*(S2z/ode.m2 - S1z/ode.m1);

    eval_R(t, f, ode, obs, Rij);
    eval_Pn(t, f, ode, obs, Pnij);
    eval_Pv(t, f, ode, obs, Pvij);    
    eval_Pnn(t, f, ode, obs, Pnnij);
    eval_Pnv(t, f, ode, obs, Pnvij);    
    eval_Pvv(t, f, ode, obs, Pvvij);    
    eval_Pnnn(t, f, ode, obs, Pnnnij);
    eval_Pnnv(t, f, ode, obs, Pnnvij);
    eval_Pnvv(t, f, ode, obs, Pnvvij);    
    eval_Pvvv(t, f, ode, obs, Pvvvij);    

    REAL hp21 = 0;
    hp21 += -2*sqrt(PI/5.0)*Rij[3*0+2];

    hp21 += ((2*sqrt(PI/5.0)*(ny*Pnij[3*0+1] - nz*Pnij[3*0+2] - nx*Pnij[3*1+1] + nx*Pnij[3*2+2] + vy*Pvij[3*0+1] - vz*Pvij[3*0+2] - vx*Pvij[3*1+1] + vx*Pvij[3*2+2]))/3.0);

    hp21 += ((-2*sqrt(PI/5.0)*(-2*nx*nz*Pnnij[3*0+0] - 6*ny*nz*Pnnij[3*0+1] + 5*nx*nx*Pnnij[3*0+2] + 11*ny*ny*Pnnij[3*0+2] + 5*nz*nz*Pnnij[3*0+2] + 4*nx*nz*Pnnij[3*1+1] - 6*nx*ny*Pnnij[3*1+2] - 2*nx*nz*Pnnij[3*2+2] - nz*vx*Pnvij[3*0+0] - nx*vz*Pnvij[3*0+0] - 3*nz*vy*Pnvij[3*0+1] - 3*ny*vz*Pnvij[3*0+1] + 5*nx*vx*Pnvij[3*0+2] + 11*ny*vy*Pnvij[3*0+2] + 5*nz*vz*Pnvij[3*0+2] + 2*nz*vx*Pnvij[3*1+1] + 2*nx*vz*Pnvij[3*1+1] - 3*ny*vx*Pnvij[3*1+2] - 3*nx*vy*Pnvij[3*1+2] - nz*vx*Pnvij[3*2+2] - nx*vz*Pnvij[3*2+2] - 2*vx*vz*Pvvij[3*0+0] - 6*vy*vz*Pvvij[3*0+1] + 5*vx*vx*Pvvij[3*0+2] + 11*vy*vy*Pvvij[3*0+2] + 5*vz*vz*Pvvij[3*0+2] + 4*vx*vz*Pvvij[3*1+1] - 6*vx*vy*Pvvij[3*1+2] - 2*vx*vz*Pvvij[3*2+2]))/21.0);

    hp21 += (-(sqrt(PI/5.0)*(3*nx*ny*ny*Pnnnij[3*0+0] - 
     3*nx*nz*nz*Pnnnij[3*0+0] - 12*nx*nx*ny*Pnnnij[3*0+1] - 
     6*ny*ny*ny*Pnnnij[3*0+1] - 12*ny*nz*nz*Pnnnij[3*0+1] + 
     12*nx*nx*nz*Pnnnij[3*0+2] + 12*ny*ny*nz*Pnnnij[3*0+2] + 
     6*nz*nz*nz*Pnnnij[3*0+2] + 9*nx*nx*nx*Pnnnij[3*1+1] + 
     6*nx*ny*ny*Pnnnij[3*1+1] + 9*nx*nz*nz*Pnnnij[3*1+1] - 
     9*nx*nx*nx*Pnnnij[3*2+2] - 9*nx*ny*ny*Pnnnij[3*2+2] - 
     6*nx*nz*nz*Pnnnij[3*2+2] + ny*ny*vx*Pnnvij[3*0+0] - 
     nz*nz*vx*Pnnvij[3*0+0] + 2*nx*ny*vy*Pnnvij[3*0+0] - 
     2*nx*nz*vz*Pnnvij[3*0+0] - 8*nx*ny*vx*Pnnvij[3*0+1] - 
     4*nx*nx*vy*Pnnvij[3*0+1] - 6*ny*ny*vy*Pnnvij[3*0+1] - 
     4*nz*nz*vy*Pnnvij[3*0+1] - 8*ny*nz*vz*Pnnvij[3*0+1] + 
     8*nx*nz*vx*Pnnvij[3*0+2] + 8*ny*nz*vy*Pnnvij[3*0+2] + 
     4*nx*nx*vz*Pnnvij[3*0+2] + 4*ny*ny*vz*Pnnvij[3*0+2] + 
     6*nz*nz*vz*Pnnvij[3*0+2] + 9*nx*nx*vx*Pnnvij[3*1+1] + 
     2*ny*ny*vx*Pnnvij[3*1+1] + 3*nz*nz*vx*Pnnvij[3*1+1] + 
     4*nx*ny*vy*Pnnvij[3*1+1] + 6*nx*nz*vz*Pnnvij[3*1+1] - 
     9*nx*nx*vx*Pnnvij[3*2+2] - 3*ny*ny*vx*Pnnvij[3*2+2] - 
     2*nz*nz*vx*Pnnvij[3*2+2] - 6*nx*ny*vy*Pnnvij[3*2+2] - 
     4*nx*nz*vz*Pnnvij[3*2+2] + 2*ny*vx*vy*Pnvvij[3*0+0] + 
     nx*vy*vy*Pnvvij[3*0+0] - 2*nz*vx*vz*Pnvvij[3*0+0] - 
     nx*vz*vz*Pnvvij[3*0+0] - 4*ny*vx*vx*Pnvvij[3*0+1] - 
     8*nx*vx*vy*Pnvvij[3*0+1] - 6*ny*vy*vy*Pnvvij[3*0+1] - 
     8*nz*vy*vz*Pnvvij[3*0+1] - 4*ny*vz*vz*Pnvvij[3*0+1] + 
     4*nz*vx*vx*Pnvvij[3*0+2] + 4*nz*vy*vy*Pnvvij[3*0+2] + 
     8*nx*vx*vz*Pnvvij[3*0+2] + 8*ny*vy*vz*Pnvvij[3*0+2] + 
     6*nz*vz*vz*Pnvvij[3*0+2] + 9*nx*vx*vx*Pnvvij[3*1+1] + 
     4*ny*vx*vy*Pnvvij[3*1+1] + 2*nx*vy*vy*Pnvvij[3*1+1] + 
     6*nz*vx*vz*Pnvvij[3*1+1] + 3*nx*vz*vz*Pnvvij[3*1+1] - 
     9*nx*vx*vx*Pnvvij[3*2+2] - 6*ny*vx*vy*Pnvvij[3*2+2] - 
     3*nx*vy*vy*Pnvvij[3*2+2] - 4*nz*vx*vz*Pnvvij[3*2+2] - 
     2*nx*vz*vz*Pnvvij[3*2+2] + 3*vx*vy*vy*Pvvvij[3*0+0] - 
     3*vx*vz*vz*Pvvvij[3*0+0] - 12*vx*vx*vy*Pvvvij[3*0+1] - 
     6*vy*vy*vy*Pvvvij[3*0+1] - 12*vy*vz*vz*Pvvvij[3*0+1] + 
     12*vx*vx*vz*Pvvvij[3*0+2] + 12*vy*vy*vz*Pvvvij[3*0+2] + 
     6*vz*vz*vz*Pvvvij[3*0+2] + 9*vx*vx*vx*Pvvvij[3*1+1] + 
     6*vx*vy*vy*Pvvvij[3*1+1] + 9*vx*vz*vz*Pvvvij[3*1+1] - 
     9*vx*vx*vx*Pvvvij[3*2+2] - 9*vx*vy*vy*Pvvvij[3*2+2] - 
     6*vx*vz*vz*Pvvvij[3*2+2]))/21.0);

    if((ode.hterms & H_PQSO) != 0) {
     hp21 += -4*ode.mu/rsq*sqrt(PI/5.0)*(ny * Dz + nz * Dy);
    }

    if((ode.hterms & H_P15QSO) != 0) {
     hp21 += 2*ode.mu*((2*sqrt(PI/5.0)*(3*dm*Dy*nx*nx*rdot - 3*dm*Dx*nx*ny*rdot + 
    3*dm*Dz*ny*nz*rdot - 3*dm*Dy*nz*nz*rdot - 3*ode.m*nx*ny*rdot*Sx + 
    3*ode.m*nx*nx*rdot*Sy - 3*ode.m*nz*nz*rdot*Sy + 3*ode.m*ny*nz*rdot*Sz - 2*dm*Dy*nx*vx - 
    2*dm*Dx*ny*vx - 2*ode.m*ny*Sx*vx - 2*ode.m*nx*Sy*vx + 4*dm*Dx*nx*vy - 
    4*dm*Dz*nz*vy + 4*ode.m*nx*Sx*vy - 4*ode.m*nz*Sz*vy + 2*dm*Dz*ny*vz + 
    2*dm*Dy*nz*vz + 2*ode.m*nz*Sy*vz + 2*ode.m*ny*Sz*vz))/(3.0*ode.m*r*r));
     hp21 += 2*ode.mu*((-4*sqrt(PI/5.0)*(dm*Dy*nx*vx - 2*dm*Dx*ny*vx - 2*ode.m*ny*Sx*vx + 
    ode.m*nx*Sy*vx + dm*Dx*nx*vy - dm*Dz*nz*vy + ode.m*nx*Sx*vy - ode.m*nz*Sz*vy + 
    2*dm*Dz*ny*vz - dm*Dy*nz*vz - ode.m*nz*Sy*vz + 2*ode.m*ny*Sz*vz))/(3.0*ode.m*r*r));
    }
   return hp21;
}

static FUNCDEF(hp2m1) {
    REAL Rij[9];
    REAL Pnij[9];
    REAL Pvij[9];
    REAL Pnnij[9];
    REAL Pnvij[9];
    REAL Pvvij[9];
    REAL Pnnnij[9];
    REAL Pnnvij[9];
    REAL Pnvvij[9];
    REAL Pvvvij[9];
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);
    REAL nx = rx/r;  
    REAL ny = ry/r;  
    REAL nz = rz/r;

    REAL vx = eval_vx(t, f, ode, obs);
    REAL vy = eval_vy(t, f, ode, obs);    
    REAL vz = eval_vz(t, f, ode, obs);  

    REAL rdot = (rx*vx + ry*vy + rz*vz)/r;
    REAL dm = ode.m1 - ode.m2;
   
    REAL cSs1 = ode.m1*ode.m1;
    REAL S1x = cSs1*f[I_s1x];
    REAL S1y = cSs1*f[I_s1y];
    REAL S1z = cSs1*f[I_s1z];
    REAL cSs2 = ode.m2*ode.m2;
    REAL S2x = cSs2*f[I_s2x];
    REAL S2y = cSs2*f[I_s2y];
    REAL S2z = cSs2*f[I_s2z];
    REAL Sx = S1x + S2x;
    REAL Sy = S1y + S2y;
    REAL Sz = S1z + S2z;
    REAL Dx = ode.m*(S2x/ode.m2 - S1x/ode.m1);
    REAL Dy = ode.m*(S2y/ode.m2 - S1y/ode.m1);
    REAL Dz = ode.m*(S2z/ode.m2 - S1z/ode.m1);

    eval_R(t, f, ode, obs, Rij);
    eval_Pn(t, f, ode, obs, Pnij);
    eval_Pv(t, f, ode, obs, Pvij);    
    eval_Pnn(t, f, ode, obs, Pnnij);
    eval_Pnv(t, f, ode, obs, Pnvij);    
    eval_Pvv(t, f, ode, obs, Pvvij);    
    eval_Pnnn(t, f, ode, obs, Pnnnij);
    eval_Pnnv(t, f, ode, obs, Pnnvij);
    eval_Pnvv(t, f, ode, obs, Pnvvij);    
    eval_Pvvv(t, f, ode, obs, Pvvvij);    

    REAL hp2m1 = 0;
    hp2m1 += 2*sqrt(PI/5.0)*Rij[3*0+2];

    hp2m1 += ((2*sqrt(PI/5.0)*(ny*Pnij[3*0+1] - nz*Pnij[3*0+2] - nx*Pnij[3*1+1] + nx*Pnij[3*2+2] + vy*Pvij[3*0+1] - vz*Pvij[3*0+2] - vx*Pvij[3*1+1] + vx*Pvij[3*2+2]))/3.0);

    hp2m1 += ((2*sqrt(PI/5.0)*(-2*nx*nz*Pnnij[3*0+0] - 6*ny*nz*Pnnij[3*0+1] + 5*nx*nx*Pnnij[3*0+2] + 11*ny*ny*Pnnij[3*0+2] + 5*nz*nz*Pnnij[3*0+2] + 4*nx*nz*Pnnij[3*1+1] - 6*nx*ny*Pnnij[3*1+2] - 2*nx*nz*Pnnij[3*2+2] - nz*vx*Pnvij[3*0+0] - nx*vz*Pnvij[3*0+0] - 3*nz*vy*Pnvij[3*0+1] - 3*ny*vz*Pnvij[3*0+1] + 5*nx*vx*Pnvij[3*0+2] + 11*ny*vy*Pnvij[3*0+2] + 5*nz*vz*Pnvij[3*0+2] + 2*nz*vx*Pnvij[3*1+1] + 2*nx*vz*Pnvij[3*1+1] - 3*ny*vx*Pnvij[3*1+2] - 3*nx*vy*Pnvij[3*1+2] - nz*vx*Pnvij[3*2+2] - nx*vz*Pnvij[3*2+2] - 2*vx*vz*Pvvij[3*0+0] - 6*vy*vz*Pvvij[3*0+1] + 5*vx*vx*Pvvij[3*0+2] + 11*vy*vy*Pvvij[3*0+2] + 5*vz*vz*Pvvij[3*0+2] + 4*vx*vz*Pvvij[3*1+1] - 6*vx*vy*Pvvij[3*1+2] - 2*vx*vz*Pvvij[3*2+2]))/21.0);

    hp2m1 += (-(sqrt(PI/5.0)*(3*nx*ny*ny*Pnnnij[3*0+0] - 
     3*nx*nz*nz*Pnnnij[3*0+0] - 12*nx*nx*ny*Pnnnij[3*0+1] - 
     6*ny*ny*ny*Pnnnij[3*0+1] - 12*ny*nz*nz*Pnnnij[3*0+1] + 
     12*nx*nx*nz*Pnnnij[3*0+2] + 12*ny*ny*nz*Pnnnij[3*0+2] + 
     6*nz*nz*nz*Pnnnij[3*0+2] + 9*nx*nx*nx*Pnnnij[3*1+1] + 
     6*nx*ny*ny*Pnnnij[3*1+1] + 9*nx*nz*nz*Pnnnij[3*1+1] - 
     9*nx*nx*nx*Pnnnij[3*2+2] - 9*nx*ny*ny*Pnnnij[3*2+2] - 
     6*nx*nz*nz*Pnnnij[3*2+2] + ny*ny*vx*Pnnvij[3*0+0] - 
     nz*nz*vx*Pnnvij[3*0+0] + 2*nx*ny*vy*Pnnvij[3*0+0] - 
     2*nx*nz*vz*Pnnvij[3*0+0] - 8*nx*ny*vx*Pnnvij[3*0+1] - 
     4*nx*nx*vy*Pnnvij[3*0+1] - 6*ny*ny*vy*Pnnvij[3*0+1] - 
     4*nz*nz*vy*Pnnvij[3*0+1] - 8*ny*nz*vz*Pnnvij[3*0+1] + 
     8*nx*nz*vx*Pnnvij[3*0+2] + 8*ny*nz*vy*Pnnvij[3*0+2] + 
     4*nx*nx*vz*Pnnvij[3*0+2] + 4*ny*ny*vz*Pnnvij[3*0+2] + 
     6*nz*nz*vz*Pnnvij[3*0+2] + 9*nx*nx*vx*Pnnvij[3*1+1] + 
     2*ny*ny*vx*Pnnvij[3*1+1] + 3*nz*nz*vx*Pnnvij[3*1+1] + 
     4*nx*ny*vy*Pnnvij[3*1+1] + 6*nx*nz*vz*Pnnvij[3*1+1] - 
     9*nx*nx*vx*Pnnvij[3*2+2] - 3*ny*ny*vx*Pnnvij[3*2+2] - 
     2*nz*nz*vx*Pnnvij[3*2+2] - 6*nx*ny*vy*Pnnvij[3*2+2] - 
     4*nx*nz*vz*Pnnvij[3*2+2] + 2*ny*vx*vy*Pnvvij[3*0+0] + 
     nx*vy*vy*Pnvvij[3*0+0] - 2*nz*vx*vz*Pnvvij[3*0+0] - 
     nx*vz*vz*Pnvvij[3*0+0] - 4*ny*vx*vx*Pnvvij[3*0+1] - 
     8*nx*vx*vy*Pnvvij[3*0+1] - 6*ny*vy*vy*Pnvvij[3*0+1] - 
     8*nz*vy*vz*Pnvvij[3*0+1] - 4*ny*vz*vz*Pnvvij[3*0+1] + 
     4*nz*vx*vx*Pnvvij[3*0+2] + 4*nz*vy*vy*Pnvvij[3*0+2] + 
     8*nx*vx*vz*Pnvvij[3*0+2] + 8*ny*vy*vz*Pnvvij[3*0+2] + 
     6*nz*vz*vz*Pnvvij[3*0+2] + 9*nx*vx*vx*Pnvvij[3*1+1] + 
     4*ny*vx*vy*Pnvvij[3*1+1] + 2*nx*vy*vy*Pnvvij[3*1+1] + 
     6*nz*vx*vz*Pnvvij[3*1+1] + 3*nx*vz*vz*Pnvvij[3*1+1] - 
     9*nx*vx*vx*Pnvvij[3*2+2] - 6*ny*vx*vy*Pnvvij[3*2+2] - 
     3*nx*vy*vy*Pnvvij[3*2+2] - 4*nz*vx*vz*Pnvvij[3*2+2] - 
     2*nx*vz*vz*Pnvvij[3*2+2] + 3*vx*vy*vy*Pvvvij[3*0+0] - 
     3*vx*vz*vz*Pvvvij[3*0+0] - 12*vx*vx*vy*Pvvvij[3*0+1] - 
     6*vy*vy*vy*Pvvvij[3*0+1] - 12*vy*vz*vz*Pvvvij[3*0+1] + 
     12*vx*vx*vz*Pvvvij[3*0+2] + 12*vy*vy*vz*Pvvvij[3*0+2] + 
     6*vz*vz*vz*Pvvvij[3*0+2] + 9*vx*vx*vx*Pvvvij[3*1+1] + 
     6*vx*vy*vy*Pvvvij[3*1+1] + 9*vx*vz*vz*Pvvvij[3*1+1] - 
     9*vx*vx*vx*Pvvvij[3*2+2] - 9*vx*vy*vy*Pvvvij[3*2+2] - 
     6*vx*vz*vz*Pvvvij[3*2+2]))/21.0);

    if((ode.hterms & H_PQSO) != 0) {
     hp2m1 += -4*ode.mu/rsq*sqrt(PI/5.0)*(ny * Dz + nz * Dy);
    }

    if((ode.hterms & H_P15QSO) != 0) {
     hp2m1 += 2*ode.mu*((-2*sqrt(PI/5.0)*(3*dm*Dy*nx*nx*rdot - 3*dm*Dx*nx*ny*rdot + 
    3*dm*Dz*ny*nz*rdot - 3*dm*Dy*nz*nz*rdot - 3*ode.m*nx*ny*rdot*Sx + 
    3*ode.m*nx*nx*rdot*Sy - 3*ode.m*nz*nz*rdot*Sy + 3*ode.m*ny*nz*rdot*Sz - 2*dm*Dy*nx*vx - 
    2*dm*Dx*ny*vx - 2*ode.m*ny*Sx*vx - 2*ode.m*nx*Sy*vx + 4*dm*Dx*nx*vy - 
    4*dm*Dz*nz*vy + 4*ode.m*nx*Sx*vy - 4*ode.m*nz*Sz*vy + 2*dm*Dz*ny*vz + 
    2*dm*Dy*nz*vz + 2*ode.m*nz*Sy*vz + 2*ode.m*ny*Sz*vz))/(3.0*ode.m*r*r));
     hp2m1 += 2*ode.mu*((4*sqrt(PI/5.0)*(dm*Dy*nx*vx - 2*dm*Dx*ny*vx - 2*ode.m*ny*Sx*vx + 
    ode.m*nx*Sy*vx + dm*Dx*nx*vy - dm*Dz*nz*vy + ode.m*nx*Sx*vy - ode.m*nz*Sz*vy + 
    2*dm*Dz*ny*vz - dm*Dy*nz*vz - ode.m*nz*Sy*vz + 2*ode.m*ny*Sz*vz))/(3.0*ode.m*r*r));
    }
   return hp2m1;
}

static FUNCDEF(hx21) {
    REAL Rij[9];
    REAL Pnij[9];
    REAL Pvij[9];
    REAL Pnnij[9];
    REAL Pnvij[9];
    REAL Pvvij[9];
    REAL Pnnnij[9];
    REAL Pnnvij[9];
    REAL Pnvvij[9];
    REAL Pvvvij[9];
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);
    REAL nx = rx/r;  
    REAL ny = ry/r;  
    REAL nz = rz/r;

    REAL vx = eval_vx(t, f, ode, obs);
    REAL vy = eval_vy(t, f, ode, obs);    
    REAL vz = eval_vz(t, f, ode, obs);  

    REAL rdot = (rx*vx + ry*vy + rz*vz)/r;
    REAL dm = ode.m1 - ode.m2;
   
    REAL cSs1 = ode.m1*ode.m1;
    REAL S1x = cSs1*f[I_s1x];
    REAL S1y = cSs1*f[I_s1y];
    REAL S1z = cSs1*f[I_s1z];
    REAL cSs2 = ode.m2*ode.m2;
    REAL S2x = cSs2*f[I_s2x];
    REAL S2y = cSs2*f[I_s2y];
    REAL S2z = cSs2*f[I_s2z];
    REAL Sx = S1x + S2x;
    REAL Sy = S1y + S2y;
    REAL Sz = S1z + S2z;
    REAL Dx = ode.m*(S2x/ode.m2 - S1x/ode.m1);
    REAL Dy = ode.m*(S2y/ode.m2 - S1y/ode.m1);
    REAL Dz = ode.m*(S2z/ode.m2 - S1z/ode.m1);

    eval_R(t, f, ode, obs, Rij);
    eval_Pn(t, f, ode, obs, Pnij);    
    eval_Pv(t, f, ode, obs, Pvij);
    eval_Pnn(t, f, ode, obs, Pnnij);    
    eval_Pnv(t, f, ode, obs, Pnvij);    
    eval_Pvv(t, f, ode, obs, Pvvij);
    eval_Pnnn(t, f, ode, obs, Pnnnij);
    eval_Pnnv(t, f, ode, obs, Pnnvij);
    eval_Pnvv(t, f, ode, obs, Pnvvij);    
    eval_Pvvv(t, f, ode, obs, Pvvvij);    

    REAL hx21 = 0;  
    hx21 += -2*sqrt(PI/5.0)*Rij[3*1+2];

    hx21 += ((-2*sqrt(PI/5.0)*(ny*Pnij[3*0+0] - nx*Pnij[3*0+1] + nz*Pnij[3*1+2] - ny*Pnij[3*2+2] + vy*Pvij[3*0+0] - vx*Pvij[3*0+1] + vz*Pvij[3*1+2] - vy*Pvij[3*2+2]))/3.0);

    hx21 += ((-2*sqrt(PI/5.0)*(4*ny*nz*Pnnij[3*0+0] - 6*nx*nz*Pnnij[3*0+1] - 6*nx*ny*Pnnij[3*0+2] - 2*ny*nz*Pnnij[3*1+1] + 11*nx*nx*Pnnij[3*1+2] + 5*ny*ny*Pnnij[3*1+2] + 5*nz*nz*Pnnij[3*1+2] - 2*ny*nz*Pnnij[3*2+2] + 2*nz*vy*Pnvij[3*0+0] + 2*ny*vz*Pnvij[3*0+0] - 3*nz*vx*Pnvij[3*0+1] - 3*nx*vz*Pnvij[3*0+1] - 3*ny*vx*Pnvij[3*0+2] - 3*nx*vy*Pnvij[3*0+2] - nz*vy*Pnvij[3*1+1] - ny*vz*Pnvij[3*1+1] + 11*nx*vx*Pnvij[3*1+2] + 5*ny*vy*Pnvij[3*1+2] + 5*nz*vz*Pnvij[3*1+2] - nz*vy*Pnvij[3*2+2] - ny*vz*Pnvij[3*2+2] + 4*vy*vz*Pvvij[3*0+0] - 6*vx*vz*Pvvij[3*0+1] - 6*vx*vy*Pvvij[3*0+2] - 2*vy*vz*Pvvij[3*1+1] + 11*vx*vx*Pvvij[3*1+2] + 5*vy*vy*Pvvij[3*1+2] + 5*vz*vz*Pvvij[3*1+2] - 2*vy*vz*Pvvij[3*2+2]))/21.0);

    hx21 += (-(sqrt(PI/5.0)*(6*nx*nx*ny*Pnnnij[3*0+0] + 
     9*ny*ny*ny*Pnnnij[3*0+0] + 9*ny*nz*nz*Pnnnij[3*0+0] - 
     6*nx*nx*nx*Pnnnij[3*0+1] - 12*nx*ny*ny*Pnnnij[3*0+1] - 
     12*nx*nz*nz*Pnnnij[3*0+1] + 3*nx*nx*ny*Pnnnij[3*1+1] - 
     3*ny*nz*nz*Pnnnij[3*1+1] + 12*nx*nx*nz*Pnnnij[3*1+2] + 
     12*ny*ny*nz*Pnnnij[3*1+2] + 6*nz*nz*nz*Pnnnij[3*1+2] - 
     9*nx*nx*ny*Pnnnij[3*2+2] - 9*ny*ny*ny*Pnnnij[3*2+2] - 
     6*ny*nz*nz*Pnnnij[3*2+2] + 4*nx*ny*vx*Pnnvij[3*0+0] + 
     2*nx*nx*vy*Pnnvij[3*0+0] + 9*ny*ny*vy*Pnnvij[3*0+0] + 
     3*nz*nz*vy*Pnnvij[3*0+0] + 6*ny*nz*vz*Pnnvij[3*0+0] - 
     6*nx*nx*vx*Pnnvij[3*0+1] - 4*ny*ny*vx*Pnnvij[3*0+1] - 
     4*nz*nz*vx*Pnnvij[3*0+1] - 8*nx*ny*vy*Pnnvij[3*0+1] - 
     8*nx*nz*vz*Pnnvij[3*0+1] + 2*nx*ny*vx*Pnnvij[3*1+1] + 
     nx*nx*vy*Pnnvij[3*1+1] - nz*nz*vy*Pnnvij[3*1+1] - 
     2*ny*nz*vz*Pnnvij[3*1+1] + 8*nx*nz*vx*Pnnvij[3*1+2] + 
     8*ny*nz*vy*Pnnvij[3*1+2] + 4*nx*nx*vz*Pnnvij[3*1+2] + 
     4*ny*ny*vz*Pnnvij[3*1+2] + 6*nz*nz*vz*Pnnvij[3*1+2] - 
     6*nx*ny*vx*Pnnvij[3*2+2] - 3*nx*nx*vy*Pnnvij[3*2+2] - 
     9*ny*ny*vy*Pnnvij[3*2+2] - 2*nz*nz*vy*Pnnvij[3*2+2] - 
     4*ny*nz*vz*Pnnvij[3*2+2] + 2*ny*vx*vx*Pnvvij[3*0+0] + 
     4*nx*vx*vy*Pnvvij[3*0+0] + 9*ny*vy*vy*Pnvvij[3*0+0] + 
     6*nz*vy*vz*Pnvvij[3*0+0] + 3*ny*vz*vz*Pnvvij[3*0+0] - 
     6*nx*vx*vx*Pnvvij[3*0+1] - 8*ny*vx*vy*Pnvvij[3*0+1] - 
     4*nx*vy*vy*Pnvvij[3*0+1] - 8*nz*vx*vz*Pnvvij[3*0+1] - 
     4*nx*vz*vz*Pnvvij[3*0+1] + ny*vx*vx*Pnvvij[3*1+1] + 
     2*nx*vx*vy*Pnvvij[3*1+1] - 2*nz*vy*vz*Pnvvij[3*1+1] - 
     ny*vz*vz*Pnvvij[3*1+1] + 4*nz*vx*vx*Pnvvij[3*1+2] + 
     4*nz*vy*vy*Pnvvij[3*1+2] + 8*nx*vx*vz*Pnvvij[3*1+2] + 
     8*ny*vy*vz*Pnvvij[3*1+2] + 6*nz*vz*vz*Pnvvij[3*1+2] - 
     3*ny*vx*vx*Pnvvij[3*2+2] - 6*nx*vx*vy*Pnvvij[3*2+2] - 
     9*ny*vy*vy*Pnvvij[3*2+2] - 4*nz*vy*vz*Pnvvij[3*2+2] - 
     2*ny*vz*vz*Pnvvij[3*2+2] + 6*vx*vx*vy*Pvvvij[3*0+0] + 
     9*vy*vy*vy*Pvvvij[3*0+0] + 9*vy*vz*vz*Pvvvij[3*0+0] - 
     6*vx*vx*vx*Pvvvij[3*0+1] - 12*vx*vy*vy*Pvvvij[3*0+1] - 
     12*vx*vz*vz*Pvvvij[3*0+1] + 3*vx*vx*vy*Pvvvij[3*1+1] - 
     3*vy*vz*vz*Pvvvij[3*1+1] + 12*vx*vx*vz*Pvvvij[3*1+2] + 
     12*vy*vy*vz*Pvvvij[3*1+2] + 6*vz*vz*vz*Pvvvij[3*1+2] - 
     9*vx*vx*vy*Pvvvij[3*2+2] - 9*vy*vy*vy*Pvvvij[3*2+2] - 
     6*vy*vz*vz*Pvvvij[3*2+2]))/21.0);

    if((ode.hterms & H_PQSO) != 0) {
     hx21 += 4*ode.mu/rsq*sqrt(PI/5.0)*(nx * Dz + nz * Dx);
    }

    if((ode.hterms & H_P15QSO) != 0) {
     hx21 += 2*ode.mu*((-2*sqrt(PI/5.0)*(-3*dm*Dy*nx*ny*rdot + 3*dm*Dx*ny*ny*rdot + 
    3*dm*Dz*nx*nz*rdot - 3*dm*Dx*nz*nz*rdot + 3*ode.m*ny*ny*rdot*Sx - 
    3*ode.m*nz*nz*rdot*Sx - 3*ode.m*nx*ny*rdot*Sy + 3*ode.m*nx*nz*rdot*Sz + 
    4*dm*Dy*ny*vx - 4*dm*Dz*nz*vx + 4*ode.m*ny*Sy*vx - 4*ode.m*nz*Sz*vx - 
    2*dm*Dy*nx*vy - 2*dm*Dx*ny*vy - 2*ode.m*ny*Sx*vy - 2*ode.m*nx*Sy*vy + 
    2*dm*Dz*nx*vz + 2*dm*Dx*nz*vz + 2*ode.m*nz*Sx*vz + 2*ode.m*nx*Sz*vz))/
  (3.0*ode.m*r*r));
     hx21 += 2*ode.mu*((4*sqrt(PI/5.0)*(dm*Dy*ny*vx - dm*Dz*nz*vx + ode.m*ny*Sy*vx - 
    ode.m*nz*Sz*vx - 2*dm*Dy*nx*vy + dm*Dx*ny*vy + ode.m*ny*Sx*vy - 2*ode.m*nx*Sy*vy + 
    2*dm*Dz*nx*vz - dm*Dx*nz*vz - ode.m*nz*Sx*vz + 2*ode.m*nx*Sz*vz))/(3.0*ode.m*r*r));
    }
    return hx21;
}

static FUNCDEF(hx2m1) {
    REAL Rij[9];
    REAL Pnij[9];
    REAL Pvij[9];
    REAL Pnnij[9];
    REAL Pnvij[9];
    REAL Pvvij[9];
    REAL Pnnnij[9];
    REAL Pnnvij[9];
    REAL Pnvvij[9];
    REAL Pvvvij[9];
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);
    REAL nx = rx/r;  
    REAL ny = ry/r;  
    REAL nz = rz/r;

    REAL vx = eval_vx(t, f, ode, obs);
    REAL vy = eval_vy(t, f, ode, obs);    
    REAL vz = eval_vz(t, f, ode, obs);  

    REAL rdot = (rx*vx + ry*vy + rz*vz)/r;
    REAL dm = ode.m1 - ode.m2;
   
    REAL cSs1 = ode.m1*ode.m1;
    REAL S1x = cSs1*f[I_s1x];
    REAL S1y = cSs1*f[I_s1y];
    REAL S1z = cSs1*f[I_s1z];
    REAL cSs2 = ode.m2*ode.m2;
    REAL S2x = cSs2*f[I_s2x];
    REAL S2y = cSs2*f[I_s2y];
    REAL S2z = cSs2*f[I_s2z];
    REAL Sx = S1x + S2x;
    REAL Sy = S1y + S2y;
    REAL Sz = S1z + S2z;
    REAL Dx = ode.m*(S2x/ode.m2 - S1x/ode.m1);
    REAL Dy = ode.m*(S2y/ode.m2 - S1y/ode.m1);
    REAL Dz = ode.m*(S2z/ode.m2 - S1z/ode.m1);

    eval_R(t, f, ode, obs, Rij);
    eval_Pn(t, f, ode, obs, Pnij);    
    eval_Pv(t, f, ode, obs, Pvij);
    eval_Pnn(t, f, ode, obs, Pnnij);    
    eval_Pnv(t, f, ode, obs, Pnvij);    
    eval_Pvv(t, f, ode, obs, Pvvij);
    eval_Pnnn(t, f, ode, obs, Pnnnij);
    eval_Pnnv(t, f, ode, obs, Pnnvij);
    eval_Pnvv(t, f, ode, obs, Pnvvij);    
    eval_Pvvv(t, f, ode, obs, Pvvvij);    

    REAL hx2m1 = 0;  
    hx2m1 += -2*sqrt(PI/5.0)*Rij[3*1+2];

    hx2m1 += ((2*sqrt(PI/5.0)*(ny*Pnij[3*0+0] - nx*Pnij[3*0+1] + nz*Pnij[3*1+2] - ny*Pnij[3*2+2] + vy*Pvij[3*0+0] - vx*Pvij[3*0+1] + vz*Pvij[3*1+2] - vy*Pvij[3*2+2]))/3.0);

    hx2m1 += ((-2*sqrt(PI/5.0)*(4*ny*nz*Pnnij[3*0+0] - 6*nx*nz*Pnnij[3*0+1] - 6*nx*ny*Pnnij[3*0+2] - 2*ny*nz*Pnnij[3*1+1] + 11*nx*nx*Pnnij[3*1+2] + 5*ny*ny*Pnnij[3*1+2] + 5*nz*nz*Pnnij[3*1+2] - 2*ny*nz*Pnnij[3*2+2] + 2*nz*vy*Pnvij[3*0+0] + 2*ny*vz*Pnvij[3*0+0] - 3*nz*vx*Pnvij[3*0+1] - 3*nx*vz*Pnvij[3*0+1] - 3*ny*vx*Pnvij[3*0+2] - 3*nx*vy*Pnvij[3*0+2] - nz*vy*Pnvij[3*1+1] - ny*vz*Pnvij[3*1+1] + 11*nx*vx*Pnvij[3*1+2] + 5*ny*vy*Pnvij[3*1+2] + 5*nz*vz*Pnvij[3*1+2] - nz*vy*Pnvij[3*2+2] - ny*vz*Pnvij[3*2+2] + 4*vy*vz*Pvvij[3*0+0] - 6*vx*vz*Pvvij[3*0+1] - 6*vx*vy*Pvvij[3*0+2] - 2*vy*vz*Pvvij[3*1+1] + 11*vx*vx*Pvvij[3*1+2] + 5*vy*vy*Pvvij[3*1+2] + 5*vz*vz*Pvvij[3*1+2] - 2*vy*vz*Pvvij[3*2+2]))/21.0);

    hx2m1 += (-(sqrt(PI/5.0)*(-6*nx*nx*ny*Pnnnij[3*0+0] - 
     9*ny*ny*ny*Pnnnij[3*0+0] - 9*ny*nz*nz*Pnnnij[3*0+0] + 
     6*nx*nx*nx*Pnnnij[3*0+1] + 12*nx*ny*ny*Pnnnij[3*0+1] + 
     12*nx*nz*nz*Pnnnij[3*0+1] - 3*nx*nx*ny*Pnnnij[3*1+1] + 
     3*ny*nz*nz*Pnnnij[3*1+1] - 12*nx*nx*nz*Pnnnij[3*1+2] - 
     12*ny*ny*nz*Pnnnij[3*1+2] - 6*nz*nz*nz*Pnnnij[3*1+2] + 
     9*nx*nx*ny*Pnnnij[3*2+2] + 9*ny*ny*ny*Pnnnij[3*2+2] + 
     6*ny*nz*nz*Pnnnij[3*2+2] - 4*nx*ny*vx*Pnnvij[3*0+0] - 
     2*nx*nx*vy*Pnnvij[3*0+0] - 9*ny*ny*vy*Pnnvij[3*0+0] - 
     3*nz*nz*vy*Pnnvij[3*0+0] - 6*ny*nz*vz*Pnnvij[3*0+0] + 
     6*nx*nx*vx*Pnnvij[3*0+1] + 4*ny*ny*vx*Pnnvij[3*0+1] + 
     4*nz*nz*vx*Pnnvij[3*0+1] + 8*nx*ny*vy*Pnnvij[3*0+1] + 
     8*nx*nz*vz*Pnnvij[3*0+1] - 2*nx*ny*vx*Pnnvij[3*1+1] - 
     nx*nx*vy*Pnnvij[3*1+1] + nz*nz*vy*Pnnvij[3*1+1] + 
     2*ny*nz*vz*Pnnvij[3*1+1] - 8*nx*nz*vx*Pnnvij[3*1+2] - 
     8*ny*nz*vy*Pnnvij[3*1+2] - 4*nx*nx*vz*Pnnvij[3*1+2] - 
     4*ny*ny*vz*Pnnvij[3*1+2] - 6*nz*nz*vz*Pnnvij[3*1+2] + 
     6*nx*ny*vx*Pnnvij[3*2+2] + 3*nx*nx*vy*Pnnvij[3*2+2] + 
     9*ny*ny*vy*Pnnvij[3*2+2] + 2*nz*nz*vy*Pnnvij[3*2+2] + 
     4*ny*nz*vz*Pnnvij[3*2+2] - 2*ny*vx*vx*Pnvvij[3*0+0] - 
     4*nx*vx*vy*Pnvvij[3*0+0] - 9*ny*vy*vy*Pnvvij[3*0+0] - 
     6*nz*vy*vz*Pnvvij[3*0+0] - 3*ny*vz*vz*Pnvvij[3*0+0] + 
     6*nx*vx*vx*Pnvvij[3*0+1] + 8*ny*vx*vy*Pnvvij[3*0+1] + 
     4*nx*vy*vy*Pnvvij[3*0+1] + 8*nz*vx*vz*Pnvvij[3*0+1] + 
     4*nx*vz*vz*Pnvvij[3*0+1] - ny*vx*vx*Pnvvij[3*1+1] - 
     2*nx*vx*vy*Pnvvij[3*1+1] + 2*nz*vy*vz*Pnvvij[3*1+1] + 
     ny*vz*vz*Pnvvij[3*1+1] - 4*nz*vx*vx*Pnvvij[3*1+2] - 
     4*nz*vy*vy*Pnvvij[3*1+2] - 8*nx*vx*vz*Pnvvij[3*1+2] - 
     8*ny*vy*vz*Pnvvij[3*1+2] - 6*nz*vz*vz*Pnvvij[3*1+2] + 
     3*ny*vx*vx*Pnvvij[3*2+2] + 6*nx*vx*vy*Pnvvij[3*2+2] + 
     9*ny*vy*vy*Pnvvij[3*2+2] + 4*nz*vy*vz*Pnvvij[3*2+2] + 
     2*ny*vz*vz*Pnvvij[3*2+2] - 6*vx*vx*vy*Pvvvij[3*0+0] - 
     9*vy*vy*vy*Pvvvij[3*0+0] - 9*vy*vz*vz*Pvvvij[3*0+0] + 
     6*vx*vx*vx*Pvvvij[3*0+1] + 12*vx*vy*vy*Pvvvij[3*0+1] + 
     12*vx*vz*vz*Pvvvij[3*0+1] - 3*vx*vx*vy*Pvvvij[3*1+1] + 
     3*vy*vz*vz*Pvvvij[3*1+1] - 12*vx*vx*vz*Pvvvij[3*1+2] - 
     12*vy*vy*vz*Pvvvij[3*1+2] - 6*vz*vz*vz*Pvvvij[3*1+2] + 
     9*vx*vx*vy*Pvvvij[3*2+2] + 9*vy*vy*vy*Pvvvij[3*2+2] + 
     6*vy*vz*vz*Pvvvij[3*2+2]))/21.0);

    if((ode.hterms & H_PQSO) != 0) {
     hx2m1 += -4*ode.mu/rsq*sqrt(PI/5.0)*(nx * Dz + nz * Dx);
    }

    if((ode.hterms & H_P15QSO) != 0) {
     hx2m1 += 2*ode.mu*((-2*sqrt(PI/5.0)*(-3*dm*Dy*nx*ny*rdot + 3*dm*Dx*ny*ny*rdot + 
    3*dm*Dz*nx*nz*rdot - 3*dm*Dx*nz*nz*rdot + 3*ode.m*ny*ny*rdot*Sx - 
    3*ode.m*nz*nz*rdot*Sx - 3*ode.m*nx*ny*rdot*Sy + 3*ode.m*nx*nz*rdot*Sz + 
    4*dm*Dy*ny*vx - 4*dm*Dz*nz*vx + 4*ode.m*ny*Sy*vx - 4*ode.m*nz*Sz*vx - 
    2*dm*Dy*nx*vy - 2*dm*Dx*ny*vy - 2*ode.m*ny*Sx*vy - 2*ode.m*nx*Sy*vy + 
    2*dm*Dz*nx*vz + 2*dm*Dx*nz*vz + 2*ode.m*nz*Sx*vz + 2*ode.m*nx*Sz*vz))/
  (3.0*ode.m*r*r));
     hx2m1 += 2*ode.mu*((4*sqrt(PI/5.0)*(dm*Dy*ny*vx - dm*Dz*nz*vx + ode.m*ny*Sy*vx - 
    ode.m*nz*Sz*vx - 2*dm*Dy*nx*vy + dm*Dx*ny*vy + ode.m*ny*Sx*vy - 2*ode.m*nx*Sy*vy + 
    2*dm*Dz*nx*vz - dm*Dx*nz*vz - ode.m*nz*Sx*vz + 2*ode.m*nx*Sz*vz))/(3.0*ode.m*r*r));
    }
    return hx2m1;
}

static FUNCDEF(hp20) {
    REAL Rij[9];
    REAL Pnnij[9];
    REAL Pnvij[9];
    REAL Pvvij[9];
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);
    REAL nx = rx/r;  
    REAL ny = ry/r;  
    REAL nz = rz/r;
   
    REAL vx = eval_vx(t, f, ode, obs);
    REAL vy = eval_vy(t, f, ode, obs);    
    REAL vz = eval_vz(t, f, ode, obs);  

    REAL rdot = (rx*vx + ry*vy + rz*vz)/r;
    REAL dm = ode.m1 - ode.m2;

    REAL cSs1 = ode.m1*ode.m1;
    REAL S1x = cSs1*f[I_s1x];
    REAL S1y = cSs1*f[I_s1y];
    REAL S1z = cSs1*f[I_s1z];
    REAL cSs2 = ode.m2*ode.m2;
    REAL S2x = cSs2*f[I_s2x];
    REAL S2y = cSs2*f[I_s2y];
    REAL S2z = cSs2*f[I_s2z];
    REAL Sx = S1x + S2x;
    REAL Sy = S1y + S2y;
    REAL Sz = S1z + S2z;
    REAL Dx = ode.m*(S2x/ode.m2 - S1x/ode.m1);
    REAL Dy = ode.m*(S2y/ode.m2 - S1y/ode.m1);
    REAL Dz = ode.m*(S2z/ode.m2 - S1z/ode.m1);

    eval_R(t, f, ode, obs, Rij);
    eval_Pnn(t, f, ode, obs, Pnnij);
    eval_Pnv(t, f, ode, obs, Pnvij);  
    eval_Pvv(t, f, ode, obs, Pvvij);

    REAL hp20 = 0;
    hp20 += -sqrt(2*PI/15.0)*(-2 * Rij[3*2+2] + Rij[3*1+1] + Rij[3*0+0]);

    hp20 += (-(sqrt(2*PI/15.0)*(nx*nx*Pnnij[3*0+0] + 5*ny*ny*Pnnij[3*0+0] + nz*nz*Pnnij[3*0+0] - 8*nx*ny*Pnnij[3*0+1] + 4*nx*nz*Pnnij[3*0+2] + 5*nx*nx*Pnnij[3*1+1] + ny*ny*Pnnij[3*1+1] + nz*nz*Pnnij[3*1+1] + 4*ny*nz*Pnnij[3*1+2] - 6*nx*nx*Pnnij[3*2+2] - 6*ny*ny*Pnnij[3*2+2] - 2*nz*nz*Pnnij[3*2+2] + nx*vx*Pnvij[3*0+0] + 5*ny*vy*Pnvij[3*0+0] + nz*vz*Pnvij[3*0+0] - 4*ny*vx*Pnvij[3*0+1] - 4*nx*vy*Pnvij[3*0+1] + 2*nz*vx*Pnvij[3*0+2] + 2*nx*vz*Pnvij[3*0+2] + 5*nx*vx*Pnvij[3*1+1] + ny*vy*Pnvij[3*1+1] + nz*vz*Pnvij[3*1+1] + 2*nz*vy*Pnvij[3*1+2] + 2*ny*vz*Pnvij[3*1+2] - 6*nx*vx*Pnvij[3*2+2] - 6*ny*vy*Pnvij[3*2+2] - 2*nz*vz*Pnvij[3*2+2] + vx*vx*Pvvij[3*0+0] + 5*vy*vy*Pvvij[3*0+0] + vz*vz*Pvvij[3*0+0] - 8*vx*vy*Pvvij[3*0+1] + 4*vx*vz*Pvvij[3*0+2] + 5*vx*vx*Pvvij[3*1+1] + vy*vy*Pvvij[3*1+1] + vz*vz*Pvvij[3*1+1] + 4*vy*vz*Pvvij[3*1+2] - 6*vx*vx*Pvvij[3*2+2] - 6*vy*vy*Pvvij[3*2+2] - 2*vz*vz*Pvvij[3*2+2]))/7.0);

    if((ode.hterms & H_P15QSO) != 0) {
     hp20 += 2*ode.mu*((-2*sqrt((2*PI)/15.0)*(3*dm*Dy*nx*nz*rdot - 3*dm*Dx*ny*nz*rdot - 
    3*ode.m*ny*nz*rdot*Sx + 3*ode.m*nx*nz*rdot*Sy - 2*dm*Dz*ny*vx - 2*dm*Dy*nz*vx - 
    2*ode.m*nz*Sy*vx - 2*ode.m*ny*Sz*vx + 2*dm*Dz*nx*vy + 2*dm*Dx*nz*vy + 
    2*ode.m*nz*Sx*vy + 2*ode.m*nx*Sz*vy))/(ode.m*r*r));
     hp20 += 2*ode.mu*((-4*sqrt((2*PI)/15.0)*(dm*Dz*ny*vx + ode.m*ny*Sz*vx - dm*Dz*nx*vy - 
    ode.m*nx*Sz*vy - dm*Dy*nx*vz + dm*Dx*ny*vz + ode.m*ny*Sx*vz - ode.m*nx*Sy*vz))/
  (ode.m*r*r));
    }
   return hp20;
}

static FUNCDEF(hx20) {
    REAL Pnij[9];
    REAL Pvij[9];
    REAL Pnnnij[9];
    REAL Pnnvij[9];
    REAL Pnvvij[9];
    REAL Pvvvij[9];
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);
    REAL nx = rx/r;  
    REAL ny = ry/r;  
    REAL nz = rz/r;
   
    REAL cSs1 = ode.m1*ode.m1;
    REAL S1x = cSs1*f[I_s1x];
    REAL S1y = cSs1*f[I_s1y];
    REAL S1z = cSs1*f[I_s1z];
    REAL cSs2 = ode.m2*ode.m2;
    REAL S2x = cSs2*f[I_s2x];
    REAL S2y = cSs2*f[I_s2y];
    REAL S2z = cSs2*f[I_s2z];
//    REAL Sx = S1x + S2x;
//    REAL Sy = S1y + S2y;
//    REAL Sz = S1z + S2z;
    REAL Dx = ode.m*(S2x/ode.m2 - S1x/ode.m1);
    REAL Dy = ode.m*(S2y/ode.m2 - S1y/ode.m1);
    REAL Dz = ode.m*(S2z/ode.m2 - S1z/ode.m1);

    REAL vx = eval_vx(t, f, ode, obs);
    REAL vy = eval_vy(t, f, ode, obs);    
    REAL vz = eval_vz(t, f, ode, obs);

//    eval_R(t, f, ode, obs, Rij);
    eval_Pn(t, f, ode, obs, Pnij);    
    eval_Pv(t, f, ode, obs, Pvij);
    eval_Pnnn(t, f, ode, obs, Pnnnij);
    eval_Pnnv(t, f, ode, obs, Pnnvij);
    eval_Pnvv(t, f, ode, obs, Pnvvij);    
    eval_Pvvv(t, f, ode, obs, Pvvvij);    

    REAL hx20 = 0;  
    hx20 += (2*sqrt((2*PI)/15.0)*(ny*Pnij[3*0+2] - nx*Pnij[3*1+2] + vy*Pvij[3*0+2] - vx*Pvij[3*1+2]));

    hx20 += (-(sqrt((2*PI)/15.0)*(3*nx*ny*nz*Pnnnij[3*0+0] - 
     3*nx*nx*nz*Pnnnij[3*0+1] + 3*ny*ny*nz*Pnnnij[3*0+1] - 
     9*nx*nx*ny*Pnnnij[3*0+2] - 9*ny*ny*ny*Pnnnij[3*0+2] - 
     6*ny*nz*nz*Pnnnij[3*0+2] - 3*nx*ny*nz*Pnnnij[3*1+1] + 
     9*nx*nx*nx*Pnnnij[3*1+2] + 9*nx*ny*ny*Pnnnij[3*1+2] + 
     6*nx*nz*nz*Pnnnij[3*1+2] + ny*nz*vx*Pnnvij[3*0+0] + 
     nx*nz*vy*Pnnvij[3*0+0] + nx*ny*vz*Pnnvij[3*0+0] - 
     2*nx*nz*vx*Pnnvij[3*0+1] + 2*ny*nz*vy*Pnnvij[3*0+1] - 
     nx*nx*vz*Pnnvij[3*0+1] + ny*ny*vz*Pnnvij[3*0+1] - 
     6*nx*ny*vx*Pnnvij[3*0+2] - 3*nx*nx*vy*Pnnvij[3*0+2] - 
     9*ny*ny*vy*Pnnvij[3*0+2] - 2*nz*nz*vy*Pnnvij[3*0+2] - 
     4*ny*nz*vz*Pnnvij[3*0+2] - ny*nz*vx*Pnnvij[3*1+1] - 
     nx*nz*vy*Pnnvij[3*1+1] - nx*ny*vz*Pnnvij[3*1+1] + 
     9*nx*nx*vx*Pnnvij[3*1+2] + 3*ny*ny*vx*Pnnvij[3*1+2] + 
     2*nz*nz*vx*Pnnvij[3*1+2] + 6*nx*ny*vy*Pnnvij[3*1+2] + 
     4*nx*nz*vz*Pnnvij[3*1+2] + nz*vx*vy*Pnvvij[3*0+0] + 
     ny*vx*vz*Pnvvij[3*0+0] + nx*vy*vz*Pnvvij[3*0+0] - 
     nz*vx*vx*Pnvvij[3*0+1] + nz*vy*vy*Pnvvij[3*0+1] - 
     2*nx*vx*vz*Pnvvij[3*0+1] + 2*ny*vy*vz*Pnvvij[3*0+1] - 
     3*ny*vx*vx*Pnvvij[3*0+2] - 6*nx*vx*vy*Pnvvij[3*0+2] - 
     9*ny*vy*vy*Pnvvij[3*0+2] - 4*nz*vy*vz*Pnvvij[3*0+2] - 
     2*ny*vz*vz*Pnvvij[3*0+2] - nz*vx*vy*Pnvvij[3*1+1] - 
     ny*vx*vz*Pnvvij[3*1+1] - nx*vy*vz*Pnvvij[3*1+1] + 
     9*nx*vx*vx*Pnvvij[3*1+2] + 6*ny*vx*vy*Pnvvij[3*1+2] + 
     3*nx*vy*vy*Pnvvij[3*1+2] + 4*nz*vx*vz*Pnvvij[3*1+2] + 
     2*nx*vz*vz*Pnvvij[3*1+2] + 3*vx*vy*vz*Pvvvij[3*0+0] - 
     3*vx*vx*vz*Pvvvij[3*0+1] + 3*vy*vy*vz*Pvvvij[3*0+1] - 
     9*vx*vx*vy*Pvvvij[3*0+2] - 9*vy*vy*vy*Pvvvij[3*0+2] - 
     6*vy*vz*vz*Pvvvij[3*0+2] - 3*vx*vy*vz*Pvvvij[3*1+1] + 
     9*vx*vx*vx*Pvvvij[3*1+2] + 9*vx*vy*vy*Pvvvij[3*1+2] + 
     6*vx*vz*vz*Pvvvij[3*1+2]))/7.0);

    if((ode.hterms & H_PQSO) != 0) {
     hx20 += 4*ode.mu/rsq*sqrt(2*PI/15.0)*(nx * Dx + ny * Dy - 2 * nz * Dz);
    }
    return hx20;
}

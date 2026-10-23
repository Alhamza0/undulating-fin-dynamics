"""
Deep analysis of MATLAB vs Python integrand difference.

The key issue: MATLAB symbolic expressions include the JACOBIAN of the
parametric surface, while Python implements f_n as the raw force 
per unit area and integrates dw*dh.

In MATLAB:
  The integration is: quad2d(f, hmin, hmax, 0, wmax)
  where w is the ACTUAL distance, not normalized.

Let's trace the MATLAB sway component f1_12 step by step.
"""
import numpy as np

def analyze_matlab_expression():
    """
    MATLAB f1_12 (sway force) = 
        (Cn*D^3*h^2*rho / 2) * sign(vn_full) * en_x_term * vn_inner^2 / delta_sq^(3/2)
    
    Python integrand =
        -0.5 * rho * Cn * en_k * vn^2 * sgn(vn)
    
    Where Python:
        vn = (D*h/delta) * vn_inner
        en_k = (D/delta) * en_x_raw
    
    So Python:
        f = -0.5 * rho * Cn * (D/delta)*en_x_raw * (D*h/delta)^2 * vn_inner^2 * sgn
          = -0.5 * rho * Cn * D^3 * h^2 * en_x_raw * vn_inner^2 / delta^3 * sgn
    
    MATLAB:
        f = (Cn*D^3*h^2*rho/2) * sign * en_x_term * vn_inner^2 / delta_sq^(3/2)
          = (Cn*D^3*h^2*rho/2) * sign * en_x_term * vn_inner^2 / delta^3
    
    So the magnitude terms match! The difference must be in how en_x_term is defined.
    
    MATLAB en_x_term uses w/wmax as w_bar:
        sin(q1) * (w/wmax - 1) - w*sin(q2)/wmax
        = sin(q1) * (w_bar - 1) - sin(q2)*w_bar
    
    Python en_x_raw:
        sin(q1) * (w_bar - 1) - sin(q2) * w_bar
    
    These are the SAME! So why the 3.27% error at w_bar=0.25?
    
    Wait — let me look more carefully at the MATLAB delta_sq expression.
    """
    
    D = 0.10
    q1, q2 = np.deg2rad(30), np.deg2rad(10)
    dq1, dq2 = 3.0, 2.5
    h = 0.02
    
    cos_dq = np.cos(q1 - q2)
    sin_dq = np.sin(q1 - q2)
    
    wmax = np.sqrt(2 * h**2 * (1 - cos_dq) + D**2)
    
    # Test at w = 0.25 * wmax
    w = 0.25 * wmax
    w_bar = w / wmax
    
    print(f"h = {h}, w = {w:.6f}, wmax = {wmax:.6f}, w_bar = {w_bar:.4f}")
    print(f"cos_dq = {cos_dq:.10f}, sin_dq = {sin_dq:.10f}")
    
    # ---- MATLAB delta_sq ----
    # MATLAB: h^2*sin_dq^2 + D^2 - D^2*w*(w/wmax - 1)*(cos_dq-1)/wmax * 2
    #
    # Let me expand: -D^2*w*(w/wmax - 1)*(cos_dq-1)/wmax * 2
    # = -D^2 * w * (w/wmax - 1) * (cos_dq - 1) * 2 / wmax
    # = -2*D^2/wmax * w * (w/wmax - 1) * (cos_dq - 1)
    # = -2*D^2 * w_bar * (w_bar - 1) * (cos_dq - 1)  (since w/wmax = w_bar)
    # = 2*D^2 * w_bar * (w_bar - 1) * (1 - cos_dq)    (flip sign)
    # 
    # Wait, (w_bar-1) is negative for w_bar < 1, and (cos_dq-1) is negative.
    # So -(negative)*(negative) = -(positive)*2 = negative? No...
    #
    # Let's compute directly:
    
    # MATLAB's delta_sq term, directly from code:
    # h^2*sin(q1-q2)^2 + D^2 - D^2*w*(w*1/wmax - 1)*(cos(q1-q2)-1)*1/wmax*2
    
    matlab_delta_sq = (h**2 * sin_dq**2 + D**2 
                       - D**2 * w * (w/wmax - 1.0) * (cos_dq - 1.0) / wmax * 2.0)
    
    # Python delta^2:
    # 2*D^2*w_bar*(w_bar-1)*(1-cos_dq) + h^2*sin_dq^2 + D^2
    python_delta_sq = 2*D**2*w_bar*(w_bar - 1)*(1 - cos_dq) + h**2*sin_dq**2 + D**2
    
    print(f"\nMATLAB delta_sq = {matlab_delta_sq:.15f}")
    print(f"Python delta_sq = {python_delta_sq:.15f}")
    print(f"Difference = {matlab_delta_sq - python_delta_sq:.2e}")
    
    # ---- MATLAB vn_inner ----
    # MATLAB: dq1 - w*(dq1*2 - cos_dq*(dq1+dq2))/wmax + (w^2*(dq1+dq2)*(cos_dq-1)) / (h^2*(cos_dq-1)*2 - D^2)
    # 
    # The denominator: h^2*(cos_dq-1)*2 - D^2 = -2*h^2*(1-cos_dq) - D^2 = -(2*h^2*(1-cos_dq) + D^2) = -wmax^2
    # So:  w^2*(dq1+dq2)*(cos_dq-1) / (-wmax^2) 
    #    = -w^2*(dq1+dq2)*(cos_dq-1) / wmax^2
    #    = w_bar^2*(dq1+dq2)*(1-cos_dq)   ... wait no
    #    = -w_bar^2 * wmax^2 * (dq1+dq2) * (cos_dq-1) / wmax^2  ... hmm
    #
    # Actually: w^2/wmax^2 = w_bar^2. And (cos_dq-1)/(-wmax^2) = (cos_dq-1)/(-wmax^2)
    # Let me just compute numerically:
    
    matlab_vn_inner = (dq1 
                       - w * (dq1*2.0 - cos_dq*(dq1+dq2)) / wmax 
                       + (w**2 * (dq1+dq2) * (cos_dq - 1.0)) / (h**2*(cos_dq-1.0)*2.0 - D**2))
    
    # Python vn_inner (before D*h/delta factor):
    python_vn_inner = (dq1
                       - (1-cos_dq)*(dq1+dq2)*w_bar**2
                       + (cos_dq*(dq1+dq2) - 2*dq1)*w_bar)
    
    print(f"\nMATLAB vn_inner = {matlab_vn_inner:.15f}")
    print(f"Python vn_inner = {python_vn_inner:.15f}")
    print(f"Difference = {matlab_vn_inner - python_vn_inner:.2e}")
    
    # ---- Expand Python vn_inner ----
    # = dq1 - (1-cos_dq)*(dq1+dq2)*w_bar^2 + (cos_dq*(dq1+dq2) - 2*dq1)*w_bar
    # 
    # MATLAB first term after dq1:
    # -w*(2*dq1 - cos_dq*(dq1+dq2))/wmax
    # = -w_bar*(2*dq1 - cos_dq*(dq1+dq2))
    # = w_bar*(cos_dq*(dq1+dq2) - 2*dq1)   ← matches Python's 3rd term
    # 
    # MATLAB second term:
    # w^2*(dq1+dq2)*(cos_dq-1) / (-wmax^2)
    # = w_bar^2 * (dq1+dq2) * (cos_dq-1) / (-1)
    # = -w_bar^2 * (dq1+dq2) * (cos_dq-1)
    # = w_bar^2 * (dq1+dq2) * (1-cos_dq)    ... wait that's wrong sign
    # 
    # Actually: denominator = h^2*(cos_dq-1)*2 - D^2 = -wmax^2
    # So: w^2*(dq1+dq2)*(cos_dq-1) / (-wmax^2)
    #   = -(w^2/wmax^2) * (dq1+dq2) * (cos_dq-1)
    #   = -w_bar^2 * (dq1+dq2) * (cos_dq-1)
    #   = w_bar^2 * (dq1+dq2) * (1-cos_dq)
    # 
    # But Python has: -(1-cos_dq)*(dq1+dq2)*w_bar^2
    # That's the NEGATIVE of what MATLAB has!
    
    print(f"\n==== DETAILED TERM BREAKDOWN ====")
    
    # MATLAB term 2 (the w*(2dq1 - cos_dq*(dq1+dq2))/wmax part):
    m_term2 = -w * (dq1*2.0 - cos_dq*(dq1+dq2)) / wmax
    # Same as: w_bar * (cos_dq*(dq1+dq2) - 2*dq1)
    p_term2 = (cos_dq*(dq1+dq2) - 2*dq1) * w_bar
    print(f"Linear w_bar term: MATLAB={m_term2:.15f}, Python={p_term2:.15f}, diff={m_term2-p_term2:.2e}")
    
    # MATLAB term 3 (the w^2 part):
    m_term3 = (w**2 * (dq1+dq2) * (cos_dq - 1.0)) / (h**2*(cos_dq-1.0)*2.0 - D**2)
    # Python equivalent: -(1-cos_dq)*(dq1+dq2)*w_bar^2
    p_term3 = -(1-cos_dq)*(dq1+dq2)*w_bar**2
    print(f"Quadratic w_bar term: MATLAB={m_term3:.15f}, Python={p_term3:.15f}, diff={m_term3-p_term3:.2e}")
    
    # AHA! The sign difference is here!
    # MATLAB: (cos_dq-1) / (h^2*(cos_dq-1)*2 - D^2) = (cos_dq-1)/(-wmax^2)
    # = -(1-cos_dq)/(-wmax^2) = (1-cos_dq)/wmax^2
    # So MATLAB term 3 = w^2*(dq1+dq2)*(1-cos_dq)/wmax^2 = +w_bar^2*(dq1+dq2)*(1-cos_dq)
    # But Python term 3 = -(1-cos_dq)*(dq1+dq2)*w_bar^2
    # => SIGN DIFFERENCE!
    
    print(f"\n❌ THE BUG IS IN THE SIGN OF THE QUADRATIC w_bar^2 TERM IN vn!")
    print(f"   MATLAB: +w_bar^2 * (dq1+dq2) * (1-cos_dq)")
    print(f"   Python: -w_bar^2 * (dq1+dq2) * (1-cos_dq)")
    print(f"   The Python sign should be POSITIVE, not negative!")
    
    # Verify with numerical expression
    correct_vn_inner = (dq1
                        + (1-cos_dq)*(dq1+dq2)*w_bar**2
                        + (cos_dq*(dq1+dq2) - 2*dq1)*w_bar)
    print(f"\n   Corrected Python vn: {correct_vn_inner:.15f}")
    print(f"   MATLAB vn:           {matlab_vn_inner:.15f}")
    print(f"   Now matches: {abs(correct_vn_inner - matlab_vn_inner) < 1e-12}")
    

if __name__ == "__main__":
    analyze_matlab_expression()

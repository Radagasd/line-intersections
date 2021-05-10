#ifndef LINE_INTERSECTIONS_INCLUDED
#define LINE_INTERSECTIONS_INCLUDED

static const float maxLineDst = 3.402823466e+38;

float3x3 constructTransitionMatrix(float3 forwardDir, float3 upDir)
{
    float3 rightDir = cross(forwardDir, upDir);
    float3x3 result = {rightDir, upDir, forwardDir};
    return result;
}

float2 quadraticFormula(float a, float b, float c)
{
    float discriminant = b * b - 4 * a * c;
    if (discriminant >= 0) 
    {
        float s = sqrt(discriminant);

        float minDst = (-b - s) / (2 * a);
        float maxDst = (-b + s) / (2 * a);
 
        return float2(minDst, maxDst);
    }
    // returns max float if no intersection
    return float2(maxLineDst, maxLineDst);
}

// Based on plane equation from https://en.wikipedia.org/wiki/Plane_(geometry)#Point%E2%80%93normal_form_and_general_form_of_the_equation_of_a_plane
float intersectPlane(float3 lineOrigin, float3 lineDir, float3 shapeOrigin, float3 shapeUpDir)
{
    // Transform line origin and direction from world space to the shape space
    float3x3 transitionMatrix = constructTransitionMatrix(float3(0,0,0), shapeUpDir);
    float3 lO = mul(transitionMatrix, lineOrigin - shapeOrigin);
    float3 lD = mul(transitionMatrix, lineDir);

    float denominator = lD.y;
    float numerator = lO.y;
    
    return - numerator / denominator;
}

// Based on quadrics equations from https://en.wikipedia.org/wiki/Quadric
float2 intersectEllipticCylinder(float3 lineOrigin, float3 lineDir, float3 shapeOrigin, float3 shapeForwardDir, float3 shapeUpDir, float2 scaleParams)
{
    // Transform line origin and direction from world space to the shape space
    float3x3 transitionMatrix = constructTransitionMatrix(shapeForwardDir, shapeUpDir);
    float3 lO = mul(transitionMatrix, lineOrigin - shapeOrigin);
    float3 lD = mul(transitionMatrix, lineDir);

    // Squared reciprocals of the shape parameters a_s and b_s
    float rsx = 1.0 / (scaleParams.x*scaleParams.x);
    float rsz = 1.0 / (scaleParams.y*scaleParams.y);

    float a = lD.x*lD.x * rsx + lD.z*lD.z * rsz;
    float b = 2 * (lO.x * lD.x * rsx + lO.z * lD.z * rsz);
    float c = lO.x*lO.x * rsx + lO.z*lO.z * rsz - 1;
    
    return quadraticFormula(a, b, c);
}

// Based on quadrics equations from https://en.wikipedia.org/wiki/Quadric
float2 intersectEllipsoid(float3 lineOrigin, float3 lineDir, float3 shapeOrigin, float3 shapeForwardDir, float3 shapeUpDir, float3 scaleParams)
{
    // Transform line origin and direction from world space to the shape space
    float3x3 transitionMatrix = constructTransitionMatrix(shapeForwardDir, shapeUpDir);
    float3 lO = mul(transitionMatrix, lineOrigin - shapeOrigin);
    float3 lD = mul(transitionMatrix, lineDir);

    // Squared reciprocals of the shape parameters
    float rsx = 1.0 / (scaleParams.x*scaleParams.x);
    float rsy = 1.0 / (scaleParams.y*scaleParams.y);
    float rsz = 1.0 / (scaleParams.z*scaleParams.z);

    float a = lD.x*lD.x * rsx + lD.y*lD.y * rsy + lD.z*lD.z * rsz;
    float b = 2 * (lO.x * lD.x * rsx + lO.y * lD.y * rsy + lO.z * lD.z * rsz);
    float c = lO.x*lO.x * rsx + lO.y*lO.y * rsy + lO.z*lO.z * rsz - 1;

    return quadraticFormula(a, b, c);
}

// Based on quadrics equations from https://en.wikipedia.org/wiki/Quadric
float2 intersectHyperboloid(float3 lineOrigin, float3 lineDir, float3 shapeOrigin, float3 shapeForwardDir, float3 shapeUpDir, float3 scaleParams, float typeParam)
{
    // Transform line origin and direction from world space to the shape space
    float3x3 transitionMatrix = constructTransitionMatrix(shapeForwardDir, shapeUpDir);
    float3 lO = mul(transitionMatrix, lineOrigin - shapeOrigin);
    float3 lD = mul(transitionMatrix, lineDir);

    // Squared reciprocals of the shape parameters
    float rsx = 1.0 / (scaleParams.x*scaleParams.x);
    float rsy = 1.0 / (scaleParams.y*scaleParams.y);
    float rsz = 1.0 / (scaleParams.z*scaleParams.z);

    float a = lD.x*lD.x * rsx - lD.y*lD.y * rsy + lD.z*lD.z * rsz;
    float b = 2 * (lO.x * lD.x * rsx - lO.y * lD.y * rsy + lO.z * lD.z * rsz);
    float c = lO.x*lO.x * rsx - lO.y*lO.y * rsy + lO.z*lO.z * rsz - typeParam;

    return quadraticFormula(a, b, c);
}

// Based on quadrics equations from https://en.wikipedia.org/wiki/Quadric
float2 intersectParaboloid(float3 lineOrigin, float3 lineDir, float3 shapeOrigin, float3 shapeForwardDir, float3 shapeUpDir, float2 scaleParams, float typeParam)
{
    // Transform line origin and direction from world space to the shape space
    float3x3 transitionMatrix = constructTransitionMatrix(shapeForwardDir, shapeUpDir);
    float3 lO = mul(transitionMatrix, lineOrigin - shapeOrigin);
    float3 lD = mul(transitionMatrix, lineDir);

    // Squared reciprocals of the shape parameters
    float rsx = 1.0 / (scaleParams.x*scaleParams.x);
    float rsz = 1.0 / (scaleParams.y*scaleParams.y);

    // Sign of type parameter
    float sst = sign(typeParam);

    float a = lD.x*lD.x * rsx + sst * lD.z*lD.z * rsz;
    float b = 2 * (lO.x * lD.x * rsx + sst * lO.z * lD.z * rsz) - lD.y ;
    float c = lO.x*lO.x * rsx + sst * lO.z*lO.z * rsz - lO.y;

    return quadraticFormula(a, b, c);
}

float intersectEllipse(float3 lineOrigin, float3 lineDir, float3 shapeOrigin, float3 shapeForwardDir, float3 shapeUpDir, float2 scaleParams)
{
    // Transform line origin and direction from world space to the shape space
    float3x3 transitionMatrix = constructTransitionMatrix(shapeForwardDir, shapeUpDir);
    float3 lO = mul(transitionMatrix, lineOrigin - shapeOrigin);
    float3 lD = mul(transitionMatrix, lineDir);

    float denominator = lD.y;
    float numerator = lO.y;

    float pIntersect = - numerator / denominator;

    float3 samplePos = lO + lD * pIntersect;
    bool isInEllipse = samplePos.x * samplePos.x / (scaleParams.x * scaleParams.x) + samplePos.z * samplePos.z / (scaleParams.y * scaleParams.y) <= 1;
    return isInEllipse ? pIntersect : maxLineDst;
}

float2 intersectCappedCylinder(float3 lineOrigin, float3 lineDir, float3 shapeOrigin, float3 shapeForwardDir, float3 shapeUpDir, float3 scaleParams)
{
    // Cylinder intersections
    float2 cylinderDst = intersectEllipticCylinder(lineOrigin, lineDir, shapeOrigin, shapeForwardDir, shapeUpDir, scaleParams.xz);

    // Cylinder intersect positions
    float3 samplePos1 = lineOrigin + lineDir * cylinderDst.x;
    float3 samplePos2 = lineOrigin + lineDir * cylinderDst.y;

    bool isMinDstBelowCap = dot(samplePos1 - shapeOrigin, shapeUpDir) <= -0.5 * scaleParams.y;
    bool isMinDstAboveCap = dot(samplePos1 - shapeOrigin, shapeUpDir) >= 0.5 * scaleParams.y;

    bool isMaxDstBelowCap = dot(samplePos2 - shapeOrigin, shapeUpDir) <= -0.5 * scaleParams.y;
    bool isMaxDstAboveCap = dot(samplePos2 - shapeOrigin, shapeUpDir) >= 0.5 * scaleParams.y;

    float minDst = cylinderDst.x;
    float maxDst = cylinderDst.y;
    if(isMinDstBelowCap || isMinDstAboveCap)
    {
        minDst = maxLineDst;
    }
    if(isMaxDstBelowCap || isMaxDstAboveCap)
    {
        maxDst = maxLineDst;
    }

    // Cap positions
    float3 ellipse1Pos = shapeOrigin - shapeUpDir * 0.5 * scaleParams.y;
    float3 ellipse2Pos = shapeOrigin + shapeUpDir * 0.5 * scaleParams.y;

    // Ellipse bottom and top cap intersections
    float ellipse1Dst = intersectEllipse(lineOrigin, lineDir, ellipse1Pos, shapeForwardDir, shapeUpDir, scaleParams.xz);
    float ellipse2Dst = intersectEllipse(lineOrigin, lineDir, ellipse2Pos, shapeForwardDir, shapeUpDir, scaleParams.xz);

    // Compare cylinder intersects to ellipse intersects for minDst
    minDst = min(min(ellipse1Dst, ellipse2Dst), minDst);
    // If value is maxLineDst make -maxLineDst for max comparison to work
    ellipse1Dst = ellipse1Dst == maxLineDst ? -maxLineDst : ellipse1Dst;
    ellipse2Dst = ellipse2Dst == maxLineDst ? -maxLineDst : ellipse2Dst;
    maxDst = maxDst == maxLineDst ? -maxLineDst : maxDst;

    // Compare cylinder intersects to ellipse intersects for maxDst
    // and make sure we don't return -maxLineDst
    maxDst = max(max(ellipse1Dst, ellipse2Dst), maxDst);
    maxDst = maxDst == -maxLineDst ? maxLineDst : maxDst;
    
    return float2(minDst, maxDst);
}

float2 intersectCappedEllipsoid(float3 lineOrigin, float3 lineDir, float3 shapeOrigin, float3 shapeForwardDir, float3 shapeUpDir, float3 scaleParams, float2 capHeight)
{
    // Ellipsoid intersections
    float2 ellipsoidDst = intersectEllipsoid(lineOrigin, lineDir, shapeOrigin, shapeForwardDir, shapeUpDir, scaleParams);

    // Cap heights in -1 to 1 range
    float cap1Height = clamp(min(capHeight.x, capHeight.y), -1, 1);
    float cap2Height = clamp(max(capHeight.x, capHeight.y), -1, 1);

    // Ellipsoid intersect positions
    float3 samplePos1 = lineOrigin + lineDir * ellipsoidDst.x;
    float3 samplePos2 = lineOrigin + lineDir * ellipsoidDst.y;

    bool isMinDstBelowCap = dot(samplePos1 - shapeOrigin, shapeUpDir) <= cap1Height * scaleParams.y;
    bool isMinDstAboveCap = dot(samplePos1 - shapeOrigin, shapeUpDir) >= cap2Height * scaleParams.y;

    bool isMaxDstBelowCap = dot(samplePos2 - shapeOrigin, shapeUpDir) <= cap1Height * scaleParams.y;
    bool isMaxDstAboveCap = dot(samplePos2 - shapeOrigin, shapeUpDir) >= cap2Height * scaleParams.y;

    float minDst = ellipsoidDst.x;
    float maxDst = ellipsoidDst.y;
    if(isMinDstBelowCap || isMinDstAboveCap)
        minDst = maxLineDst;
    if(isMaxDstBelowCap || isMaxDstAboveCap)
        maxDst = maxLineDst;
    
    // Cap positions
    float3 ellipse1Pos = shapeOrigin + cap1Height * scaleParams.y * shapeUpDir;
    float3 ellipse2Pos = shapeOrigin + cap2Height * scaleParams.y * shapeUpDir;

    // Ellipsoidal cap radius from https://keisan.casio.com/keisan/image/volume%20of%20an%20ellipsoidal%20cap.pdf
    float e1HeightSqr = cap1Height * cap1Height;
    float e1ScaleX = scaleParams.x * sqrt(1 - e1HeightSqr);
    float e1ScaleZ = scaleParams.z * sqrt(1 - e1HeightSqr);

    float e2HeightSqr = cap2Height * cap2Height;
    float e2ScaleX = scaleParams.x * sqrt(1 - e2HeightSqr);
    float e2ScaleZ = scaleParams.z * sqrt(1 - e2HeightSqr);

    // Ellipse bottom and top cap intersections
    float ellipse1Dst = intersectEllipse(lineOrigin, lineDir, ellipse1Pos, shapeForwardDir, shapeUpDir, float2(e1ScaleX, e1ScaleZ));
    float ellipse2Dst = intersectEllipse(lineOrigin, lineDir, ellipse2Pos, shapeForwardDir, shapeUpDir, float2(e2ScaleX, e2ScaleZ));

    // Compare ellipsoid intersects to ellipse intersects for minDst
    minDst = min(min(ellipse1Dst, ellipse2Dst), minDst);
    // If value is maxLineDst make -maxLineDst for max comparison to work
    ellipse1Dst = ellipse1Dst == maxLineDst ? -maxLineDst : ellipse1Dst;
    ellipse2Dst = ellipse2Dst == maxLineDst ? -maxLineDst : ellipse2Dst;
    maxDst = maxDst == maxLineDst ? -maxLineDst : maxDst;

    // Compare ellipsoid intersects to ellipse intersects for maxDst
    // and make sure we don't return -maxLineDst
    maxDst = max(max(ellipse1Dst, ellipse2Dst), maxDst);
    maxDst = maxDst == -maxLineDst ? maxLineDst : maxDst;
    
    return float2(minDst,maxDst);
}

float2 intersectCappedHalfHyperboloid(float3 lineOrigin, float3 lineDir, float3 shapeOrigin, float3 shapeForwardDir, float3 shapeUpDir, float3 scaleParams, float typeParam, float2 capHeight)
{
    // Hyperboloid intersections
    float2 hyperboloidDst = intersectHyperboloid(lineOrigin, lineDir, shapeOrigin, shapeForwardDir, shapeUpDir, scaleParams, typeParam);

    // Cap heights in -maxLineDst to 0 range
    float cap1Height = clamp(min(capHeight.x, capHeight.y), -maxLineDst, 0);
    float cap2Height = clamp(max(capHeight.x, capHeight.y), -maxLineDst, 0);

    // Ellipsoid intersect positions
    float3 samplePos1 = lineOrigin + lineDir * hyperboloidDst.x;
    float3 samplePos2 = lineOrigin + lineDir * hyperboloidDst.y;

    bool isMinDstBelowCap = dot(samplePos1 - shapeOrigin, shapeUpDir) <= cap1Height;
    bool isMinDstAboveCap = dot(samplePos1 - shapeOrigin, shapeUpDir) >= cap2Height;

    bool isMaxDstBelowCap = dot(samplePos2 - shapeOrigin, shapeUpDir) <= cap1Height;
    bool isMaxDstAboveCap = dot(samplePos2 - shapeOrigin, shapeUpDir) >= cap2Height;

    float minDst = hyperboloidDst.x;
    float maxDst = hyperboloidDst.y;
    if(isMinDstBelowCap || isMinDstAboveCap)
        minDst = maxLineDst;
    if(isMaxDstBelowCap || isMaxDstAboveCap)
        maxDst = maxLineDst;
    
    // Cap positions
    float3 ellipse1Pos = shapeOrigin + cap1Height * shapeUpDir;
    float3 ellipse2Pos = shapeOrigin + cap2Height * shapeUpDir;

    // Hyperboloid cap radius
    float sySqrR = 1 / (scaleParams.y * scaleParams.y);
    float e1HeightSqr = cap1Height * cap1Height;
    float e1ScaleX = scaleParams.x * sqrt(typeParam + e1HeightSqr * sySqrR);
    float e1ScaleZ = scaleParams.z * sqrt(typeParam + e1HeightSqr * sySqrR);

    float e2HeightSqr = cap2Height * cap2Height;
    float e2ScaleX = scaleParams.x * sqrt(typeParam + e2HeightSqr * sySqrR);
    float e2ScaleZ = scaleParams.z * sqrt(typeParam + e2HeightSqr * sySqrR);

    // Ellipse bottom and top cap intersections
    float ellipse1Dst = intersectEllipse(lineOrigin, lineDir, ellipse1Pos, shapeForwardDir, shapeUpDir, float2(e1ScaleX, e1ScaleZ));
    float ellipse2Dst = intersectEllipse(lineOrigin, lineDir, ellipse2Pos, shapeForwardDir, shapeUpDir, float2(e2ScaleX, e2ScaleZ));

    // Compare ellipsoid intersects to ellipse intersects for minDst
    minDst = min(min(ellipse1Dst, ellipse2Dst), minDst);
    // If value is maxLineDst make -maxLineDst for max comparison to work
    ellipse1Dst = ellipse1Dst == maxLineDst ? -maxLineDst : ellipse1Dst;
    ellipse2Dst = ellipse2Dst == maxLineDst ? -maxLineDst : ellipse2Dst;
    maxDst = maxDst == maxLineDst ? -maxLineDst : maxDst;

    // Compare ellipsoid intersects to ellipse intersects for maxDst
    // and make sure we don't return -maxLineDst
    maxDst = max(max(ellipse1Dst, ellipse2Dst), maxDst);
    maxDst = maxDst == -maxLineDst ? maxLineDst : maxDst;
    
    return float2(minDst,maxDst);
}

float intersectRectangle(float3 lineOrigin, float3 lineDir, float3 shapeOrigin, float3 shapeForwardDir, float3 shapeUpDir, float2 scaleParams)
{
    // Transform line origin and direction from world space to the shape space
    float3x3 transitionMatrix = constructTransitionMatrix(shapeForwardDir, shapeUpDir);
    float3 lO = mul(transitionMatrix, lineOrigin - shapeOrigin);
    float3 lD = mul(transitionMatrix, lineDir);

    float denominator = lD.y;
    float numerator = lO.y;

    float pIntersect = - numerator / denominator;

    float3 samplePos = lO + lD * pIntersect;
    bool isInRectangle = abs(samplePos.x) <= scaleParams.x && abs(samplePos.z) <= scaleParams.y;
    return isInRectangle ? pIntersect : maxLineDst;
}

float2 intersectBox(float3 lineOrigin, float3 lineDir, float3 shapeOrigin, float3 shapeForwardDir, float3 shapeUpDir, float3 scaleParams)
{
    // Calculate right direction
    float3 shapeRightDir = cross(shapeForwardDir, shapeUpDir);

    // Rectangle face intersections
    float faceUp = intersectRectangle(lineOrigin, lineDir, shapeOrigin + shapeUpDir * scaleParams.y, shapeForwardDir, shapeUpDir, scaleParams.xz);
    float faceDown = intersectRectangle(lineOrigin, lineDir, shapeOrigin - shapeUpDir * scaleParams.y, shapeForwardDir, shapeUpDir, scaleParams.xz);

    float faceForward = intersectRectangle(lineOrigin, lineDir, shapeOrigin + shapeForwardDir * scaleParams.z, -shapeUpDir, shapeForwardDir, scaleParams.xy);
    float faceBackward = intersectRectangle(lineOrigin, lineDir, shapeOrigin - shapeForwardDir * scaleParams.z, -shapeUpDir, shapeForwardDir, scaleParams.xy);

    float faceRight = intersectRectangle(lineOrigin, lineDir, shapeOrigin + shapeRightDir * scaleParams.x, shapeForwardDir, shapeRightDir, scaleParams.yz);
    float faceLeft = intersectRectangle(lineOrigin, lineDir, shapeOrigin - shapeRightDir * scaleParams.x, shapeForwardDir, shapeRightDir, scaleParams.yz);

    // Min face comparison
    float minDst = min(faceUp,min(faceDown,min(faceForward,min(faceBackward, min(faceRight, faceLeft)))));
    
    // If value is maxLineDst make -maxLineDst for max comparison to work
    faceUp *= faceUp == maxLineDst ? -1 : 1;
    faceDown *= faceDown == maxLineDst ? -1 : 1;
    faceForward *= faceForward == maxLineDst ? -1 : 1;
    faceBackward *= faceBackward == maxLineDst ? -1 : 1;
    faceRight *= faceRight == maxLineDst ? -1 : 1;
    faceLeft *= faceLeft == maxLineDst ? -1 : 1;

    // Maxn face comparison and make sure we don't return -maxLineDst
    float maxDst = max(faceUp,max(faceDown,max(faceForward,max(faceBackward, max(faceRight, faceLeft)))));
    maxDst *= maxDst == -maxLineDst ? -1 : 1;

    return float2(minDst, maxDst);
}
#endif
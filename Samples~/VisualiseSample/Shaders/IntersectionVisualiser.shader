Shader "KelvinvanHoorn/IntersectionVisualiser"
{
    Properties
    {
        _ShapeParams ("Shape parameters", vector) = (1,1,1,1)
        _DistanceScale ("Distance scale", float) = 100
        _CapHeights ("Cap heights", vector) = (0,1,0,0)
        [KeywordEnum(Plane, Cylinder, Ellipsoid, Hyperboloid, Paraboloid, Ellipse, CappedCylinder, CappedEllipsoid, CappedHalfHyperboloid, Rectangle, Box)] _Shape ("Intersect shape", float) = 0
    }
    SubShader
    {
        Tags { "RenderType" = "Opaque" "RenderPipeline" = "UniversalRenderPipeline" "Queue" = "Geometry" }
        Cull Front
 
        Pass
        {
            HLSLPROGRAM
            #pragma vertex vert
            #pragma fragment frag

            #pragma multi_compile _SHAPE_PLANE _SHAPE_CYLINDER _SHAPE_ELLIPSOID _SHAPE_HYPERBOLOID _SHAPE_PARABOLOID _SHAPE_ELLIPSE _SHAPE_CAPPEDCYLINDER _SHAPE_CAPPEDELLIPSOID _SHAPE_CAPPEDHALFHYPERBOLOID _SHAPE_RECTANGLE _SHAPE_BOX
 
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"    
            #include "Packages/com.kelvinvanhoorn.line-intersections/LineIntersections.hlsl"   
 
            struct Attributes
            {
                float4 posOS    : POSITION;
            };
 
            struct v2f
            {
                float4 posCS        : SV_POSITION;
                float3 posWS        : TEXCOORD0;
 
                float3 objectOrigin       : TEXCOORD1;
                float3 objectScale  : TEXCOORD2;
            };

            float4 _ShapeParams;
            float _DistanceScale;
            float2 _CapHeights;
 
            v2f vert(Attributes IN)
            {
                v2f OUT = (v2f)0;
 
                VertexPositionInputs vertexInput = GetVertexPositionInputs(IN.posOS.xyz);
 
                OUT.posCS = vertexInput.positionCS;
                OUT.posWS = vertexInput.positionWS;
 
                // Object information
                OUT.objectOrigin = UNITY_MATRIX_M._m03_m13_m23;
 
                return OUT;
            }
 
            float4 frag (v2f IN) : SV_Target
            {
                // Initial ray information
                float3 lineOrigin = _WorldSpaceCameraPos;
                float3 lineDir = normalize(IN.posWS - _WorldSpaceCameraPos);

                // Shape information
                float3 shapeUpDir = normalize(mul(unity_ObjectToWorld, float4(0,1,0,0)).xyz);
                float3 shapeForwardDir = normalize(mul(unity_ObjectToWorld, float4(0,0,1,0)).xyz);
                float3 shapeOrigin = IN.objectOrigin;

                float r = 0;
                float g = 0;
                float b = 0;

                #ifdef _SHAPE_PLANE
                    float intersect = intersectPlane(lineOrigin, lineDir, shapeOrigin, shapeUpDir);
                    if(intersect < maxLineDst && intersect > 0)
                    {
                        float3 pos = lineOrigin - shapeOrigin + lineDir * intersect;
                        r = saturate(abs(pos.x)/_DistanceScale);
                        g = saturate(abs(pos.y)/_DistanceScale);
                        b = saturate(abs(pos.z)/_DistanceScale);
                    }
                #endif
                #ifdef _SHAPE_CYLINDER
                    float2 intersect = intersectEllipticCylinder(lineOrigin, lineDir, shapeOrigin, shapeForwardDir, shapeUpDir, _ShapeParams.xz);
                    float frontIntersection = intersect.x >= 0 ? intersect.x : intersect.y;
                    if(frontIntersection < maxLineDst && frontIntersection >= 0)
                    {
                        float3 pos = lineOrigin - shapeOrigin + lineDir * frontIntersection;
                        r = saturate(abs(pos.x)/_DistanceScale);
                        g = saturate(abs(pos.y)/_DistanceScale);
                        b = saturate(abs(pos.z)/_DistanceScale);
                    }
                #endif
                #ifdef _SHAPE_ELLIPSOID
                    float2 intersect = intersectEllipsoid(lineOrigin, lineDir, shapeOrigin, shapeForwardDir, shapeUpDir, _ShapeParams.xyz);
                    float frontIntersection = intersect.x >= 0 ? intersect.x : intersect.y;
                    if(frontIntersection < maxLineDst && frontIntersection >= 0)
                    {
                        float3 pos = lineOrigin - shapeOrigin + lineDir * frontIntersection;
                        r = saturate(abs(pos.x)/_DistanceScale);
                        g = saturate(abs(pos.y)/_DistanceScale);
                        b = saturate(abs(pos.z)/_DistanceScale);
                    }
                #endif
                #ifdef _SHAPE_HYPERBOLOID
                    float2 intersect = intersectHyperboloid(lineOrigin, lineDir, shapeOrigin, shapeForwardDir, shapeUpDir, _ShapeParams.xyz, _ShapeParams.w);
                    float frontIntersection = intersect.x >= 0 ? intersect.x : intersect.y;
                    if(frontIntersection < maxLineDst && frontIntersection >= 0)
                    {
                        float3 pos = lineOrigin - shapeOrigin + lineDir * frontIntersection;
                        r = saturate(abs(pos.x)/_DistanceScale);
                        g = saturate(abs(pos.y)/_DistanceScale);
                        b = saturate(abs(pos.z)/_DistanceScale);
                    }
                #endif
                #ifdef _SHAPE_PARABOLOID
                    float2 intersect = intersectParaboloid(lineOrigin, lineDir, shapeOrigin, shapeForwardDir, shapeUpDir, _ShapeParams.xz, _ShapeParams.w);
                    float frontIntersection = intersect.x >= 0 ? intersect.x : intersect.y;
                    if(frontIntersection < maxLineDst && frontIntersection >= 0)
                    {
                        float3 pos = lineOrigin - shapeOrigin + lineDir * frontIntersection;
                        r = saturate(abs(pos.x)/_DistanceScale);
                        g = saturate(abs(pos.y)/_DistanceScale);
                        b = saturate(abs(pos.z)/_DistanceScale);
                    }
                #endif
                #ifdef _SHAPE_ELLIPSE
                    float intersect = intersectEllipse(lineOrigin, lineDir, shapeOrigin, shapeForwardDir, shapeUpDir, _ShapeParams.xz);
                    if(intersect < maxLineDst && intersect > 0)
                    {
                        float3 pos = lineOrigin - shapeOrigin + lineDir * intersect;
                        r = saturate(abs(pos.x)/_DistanceScale);
                        g = saturate(abs(pos.y)/_DistanceScale);
                        b = saturate(abs(pos.z)/_DistanceScale);
                    }
                #endif
                #ifdef _SHAPE_CAPPEDCYLINDER
                    float2 intersect = intersectCappedCylinder(lineOrigin, lineDir, shapeOrigin, shapeForwardDir, shapeUpDir, _ShapeParams.xyz);
                    float frontIntersection = intersect.x >= 0 ? intersect.x : intersect.y;
                    if(frontIntersection < maxLineDst && frontIntersection >= 0)
                    {
                        float3 pos = lineOrigin - shapeOrigin + lineDir * frontIntersection;
                        r = saturate(abs(pos.x)/_DistanceScale);
                        g = saturate(abs(pos.y)/_DistanceScale);
                        b = saturate(abs(pos.z)/_DistanceScale);
                    }
                #endif
                #ifdef _SHAPE_CAPPEDELLIPSOID
                    float2 intersect = intersectCappedEllipsoid(lineOrigin, lineDir, shapeOrigin, shapeForwardDir, shapeUpDir, _ShapeParams.xyz, _CapHeights);
                    float frontIntersection = intersect.x >= 0 ? intersect.x : intersect.y;
                    if(frontIntersection < maxLineDst && frontIntersection >= 0)
                    {
                        float3 pos = lineOrigin - shapeOrigin + lineDir * frontIntersection;
                        r = saturate(abs(pos.x)/_DistanceScale);
                        g = saturate(abs(pos.y)/_DistanceScale);
                        b = saturate(abs(pos.z)/_DistanceScale);
                    }
                #endif
                #ifdef _SHAPE_CAPPEDHALFHYPERBOLOID
                    float2 intersect = intersectCappedHalfHyperboloid(lineOrigin, lineDir, shapeOrigin, shapeForwardDir, shapeUpDir, _ShapeParams.xyz, _ShapeParams.w, _CapHeights);
                    float frontIntersection = intersect.x >= 0 ? intersect.x : intersect.y;
                    if(frontIntersection < maxLineDst && frontIntersection >= 0)
                    {
                        float3 pos = lineOrigin - shapeOrigin + lineDir * frontIntersection;
                        r = saturate(abs(pos.x)/_DistanceScale);
                        g = saturate(abs(pos.y)/_DistanceScale);
                        b = saturate(abs(pos.z)/_DistanceScale);
                    }
                #endif    
                #ifdef _SHAPE_RECTANGLE
                    float intersect = intersectRectangle(lineOrigin, lineDir, shapeOrigin, shapeForwardDir, shapeUpDir, _ShapeParams.xz);
                    if(intersect < maxLineDst && intersect > 0)
                    {
                        float3 pos = lineOrigin - shapeOrigin + lineDir * intersect;
                        r = saturate(abs(pos.x)/_DistanceScale);
                        g = saturate(abs(pos.y)/_DistanceScale);
                        b = saturate(abs(pos.z)/_DistanceScale);
                    }
                #endif
                #ifdef _SHAPE_BOX
                    float2 intersect = intersectBox(lineOrigin, lineDir, shapeOrigin, shapeForwardDir, shapeUpDir, _ShapeParams.xyz);
                    float frontIntersection = intersect.x >= 0 ? intersect.x : intersect.y;
                    if(frontIntersection < maxLineDst && frontIntersection >= 0)
                    {
                        float3 pos = lineOrigin - shapeOrigin + lineDir * frontIntersection;
                        r = saturate(abs(pos.x)/_DistanceScale);
                        g = saturate(abs(pos.y)/_DistanceScale);
                        b = saturate(abs(pos.z)/_DistanceScale);
                    }
                #endif
     
                return float4(r, g, b, 1);
            }
            ENDHLSL
        }
    }
}
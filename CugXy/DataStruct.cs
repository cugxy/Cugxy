using System;
using LinearAlgebra;
using System.Drawing;
using System.Collections.Generic;
using ESRI.ArcGIS.Geodatabase;
using ESRI.ArcGIS.Geometry;

namespace CugXy
{
    //struct MyLine
    //{
    //    public List<PointF> points;
    //    public double length;
    //    public MyLine(List<PointF> points, double len)
    //    {
    //        this.points = points;
    //        this.length = len;
    //    }
    //}
    struct MyPoint
    {
        public IPoint point;
        public List<IPolyline> lines;
        public MyPoint(IPoint p, List<IPolyline> lines)
        {
            this.point = p;
            this.lines = lines;
        }
    }
    class DataStruct
    {
        public List<MyPoint> MyPoints;
        public DataStruct(IGeometricNetwork pGeomNet)
        {
            MyPoints = new List<MyPoint>();
            if (pGeomNet == null)
            {
                Console.WriteLine("IGeometricNetwork error datastruct init fail");
                return;
            }
            INetwork pNet = pGeomNet.Network;
            if (pNet == null)
            {
                Console.WriteLine("network error datastruct init fail");
                return;
            }
            int edgeCount = pNet.EdgeCount;
            int junctionCount = pNet.JunctionCount;
            List<IPoint> points = new List<IPoint>(); //点集
            for (int i = 1; i <= junctionCount; i++)
            {
                IGeometry pGeom = pGeomNet.GeometryForJunctionEID[i];
                if (pGeom != null && pGeom.GeometryType == esriGeometryType.esriGeometryPoint)
                {
                    points.Add(pGeom as IPoint);
                }
            }
            List<IPolyline> lines = new List<IPolyline>(); //获取所有边
            for (int j = 1; j <= edgeCount; j++)
            {
                IGeometry pGeom = pGeomNet.GeometryForEdgeEID[j];
                if (pGeom != null && pGeom.GeometryType == esriGeometryType.esriGeometryPolyline)
                {
                    lines.Add(pGeom as IPolyline);
                }
            }
            foreach (var onePoint in points)
            {
                List<IPolyline> connEdges = new List<IPolyline>();
                foreach (var oneEdge in lines)
                {
                    IPoint fromPoint = oneEdge.FromPoint;
                    IPoint toPoint = oneEdge.ToPoint;

                    if ((fromPoint.X == onePoint.X && fromPoint.Y == onePoint.Y) || (toPoint.X == onePoint.X && toPoint.Y == onePoint.Y))
                    {
                        connEdges.Add(oneEdge);
                    }
                }
                MyPoint myPoint = new MyPoint(onePoint, connEdges);
                MyPoints.Add(myPoint);
            }
        }
    }
}


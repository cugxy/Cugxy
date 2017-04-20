using System;
using System.Collections.Generic;
using System.IO;
using System.Drawing;
using VoronoiCSharpLib;
using LinearAlgebra;
//using Esri.FileGDB;
//using OSGeo.GDAL;
//using OSGeo.OGR;
using ESRI.ArcGIS.Geodatabase;
using ESRI.ArcGIS.Controls;
using ESRI.ArcGIS.Geometry;
using ESRI.ArcGIS.DataSourcesGDB;
using ESRI.ArcGIS.NetworkAnalysis;
using ESRI.ArcGIS.Carto;
using ESRI.ArcGIS.Display;
using ESRI.ArcGIS.esriSystem;
using ESRI.ArcGIS;

namespace CSharpTest
{

    class Program
    {
        ////测试 gdal for c#
        //static void test()
        //{
        //    Ogr.RegisterAll();
        //    DataSource data = Ogr.Open("E:\\data\\test_1.gdb", 0);
        //    Layer layer = data.GetLayerByName("first1");
        //    long geomCount = layer.GetFeatureCount(0);
        //    for (int i = 1; i <= geomCount; i++)
        //    {
        //        OSGeo.OGR.Feature feature = layer.GetFeature(i);
        //        int ID = feature.GetFieldAsInteger(6);
        //        int enumId = feature.GetFieldAsInteger(7);
        //        Geometry geom = feature.GetGeometryRef();
        //        string name = geom.GetGeometryName();
        //        Geometry g = geom.GetCurveGeometry(null);
        //        int a = geom.GetGeometryCount();
        //        if (geom.GetGeometryType() == OSGeo.OGR.wkbGeometryType.wkbMultiLineString)
        //        {
        //            int pointCount = geom.GetPointCount();
        //            double length = geom.Length();

        //        }
        //    }
        //}
        ////读取文件,获取结点信息by fileGDBAPI
        //static List<Esri.FileGDB.Point> ReadFileToGetPoint(String filePath, String tablePath)
        //{
        //    try
        //    {
        //        List<Esri.FileGDB.Point> points = new List<Esri.FileGDB.Point>();
        //        Geodatabase geodatabase = Geodatabase.Open(filePath);
        //        Esri.FileGDB.Table netTable = geodatabase.OpenTable(tablePath);
        //        String ruleStr = "";
        //        RowCollection attrQueryRows = netTable.Search("SHAPE", ruleStr, RowInstance.Recycle);
        //        foreach (var oneRow in attrQueryRows)
        //        {
        //            PointShapeBuffer geometry = oneRow.GetGeometry();
        //            Esri.FileGDB.Point point = geometry.point;
        //            points.Add(point);
        //        }
        //        netTable.Close();
        //        geodatabase.Close();
        //        return points;
        //    }
        //    catch (FileGDBException ex)
        //    {
        //        Console.WriteLine("{0} - {1}", ex.Message, ex.ErrorCode);
        //        return null;
        //    }
        //    catch (Exception ex)
        //    {
        //        Console.WriteLine("General exception.  " + ex.Message);
        //        return null;
        //    }
        //}

        ////读取文件，获取结点信息By gdal
        //static Dictionary<int, PointF> GetPoint(string filePath, string layerName)
        //{
        //    Dictionary<int, PointF> result = new Dictionary<int, PointF>();
        //    Ogr.RegisterAll();
        //    DataSource dataSouce = Ogr.Open(filePath, 0);
        //    if (dataSouce == null)
        //    {
        //        Console.WriteLine("data souce error");
        //        return null;
        //    }
        //    Layer layer = dataSouce.GetLayerByName(layerName);
        //    if (layer == null)
        //    {
        //        Console.WriteLine("layer error");
        //        return null;
        //    }
        //    long count = layer.GetFeatureCount(0);
        //    int id = 0;
        //    for (int i = 1; i <= count; i++)
        //    {
        //        OSGeo.OGR.Feature feature = layer.GetFeature(i);
        //        if (feature == null)
        //        {
        //            Console.WriteLine("feature error.  i =" + i);
        //            continue;
        //        }
        //        //int Id = feature.GetFieldAsInteger(1);
        //        Geometry geom = feature.GetGeometryRef();
        //        if (geom == null)
        //        {
        //            Console.WriteLine("geometry error.  i =" + i);
        //            continue;
        //        }
        //        if (geom.GetGeometryType() == OSGeo.OGR.wkbGeometryType.wkbPoint)
        //        {
        //            double[] point = new double[3];
        //            geom.GetPoint(0, point);
        //            PointF pointf = new PointF((float)point[0], (float)point[1]);
        //            result.Add(id, pointf);
        //            id++;
        //        }
        //    }
        //    return result;
        //}

        //读取几何网络，获取结点及边的信息 by arcgis engine
        static IGeometricNetwork ReadNet(string filePath, string dataSetName, string netName)
        {
            IWorkspaceFactory pWF = new FileGDBWorkspaceFactoryClass();
            IFeatureWorkspace space = (IFeatureWorkspace)pWF.OpenFromFile(filePath, 0);
            if (space == null)
            {
                Console.WriteLine("featrue workspace error");
                return null;
            }
            IFeatureDataset dataset = space.OpenFeatureDataset(dataSetName);
            if (dataset == null)
            {
                Console.WriteLine("feature dataset error");
                return null;
            }
            INetworkCollection pNetColl = dataset as INetworkCollection;
            IGeometricNetwork pGeomNet = pNetColl.GeometricNetworkByName[netName];
            if (pGeomNet == null)
            {
                Console.WriteLine("geometry network error");
                return null;
            }
            return pGeomNet;
            // INetwork network = pGeomNet.Network;
            // int edgeCount = network.EdgeCount;
            // int junctionCount = network.JunctionCount;

            // IGeometry pEdge =  pGeomNet.GeometryForEdgeEID[1];
            // IGeometry pPoint = pGeomNet.GeometryForJunctionEID[1];

            //if(pPoint.GeometryType == esriGeometryType.esriGeometryPoint)
            // {
            //     IPoint ipoint = (pPoint as IPoint);
            // }
        }

        static List<Dictionary<int, PointF>> OneToOneMatch(IGeometricNetwork pGeomNet1, IGeometricNetwork pGeomNet2, double deviation = 3000.0)
        {
            List<Dictionary<int, PointF>> result = new List<Dictionary<int, PointF>>();
            Dictionary<int, PointF> result1 = new Dictionary<int, PointF>();
            Dictionary<int, PointF> result2 = new Dictionary<int, PointF>();
            INetwork pNet1 = pGeomNet1.Network;
            if (pNet1 == null)
            {
                Console.WriteLine("network error");
                return null;
            }
            int edgeCount1 = pNet1.EdgeCount;
            int junctionCount1 = pNet1.JunctionCount;
            List<IPoint> points1 = new List<IPoint>();
            for (int i1 = 1; i1 <= junctionCount1; i1++)
            {
                IGeometry pGeom = pGeomNet1.GeometryForJunctionEID[i1];
                if (pGeom != null && pGeom.GeometryType == esriGeometryType.esriGeometryPoint)
                {
                    points1.Add(pGeom as IPoint);
                }
            }

            INetwork pNet2 = pGeomNet2.Network;
            if (pNet1 == null)
            {
                Console.WriteLine("network error");
                return null;
            }
            int edgeCount2 = pNet2.EdgeCount;
            int junctionCount2 = pNet2.JunctionCount;
            List<IPoint> points2 = new List<IPoint>();
            for (int i2 = 1; i2 <= junctionCount2; i2++)
            {
                IGeometry pGeom = pGeomNet2.GeometryForJunctionEID[i2];
                if (pGeom != null && pGeom.GeometryType == esriGeometryType.esriGeometryPoint)
                {
                    points2.Add(pGeom as IPoint);
                }
            }

            if (junctionCount1 < junctionCount2)
            {
                int id = 0;
                foreach (var onePointInNet1 in points1)
                {
                    //pGeomNet2.JunctionElement[onePoint];
                    //如何获取该点位置在第二张网络中附近的点呢
                    //如果在第二张网络中遍历所有的点 判断距离 会不会太慢
                    List<IPoint> mayMatchPointsInNet2 = new List<IPoint>();
                    foreach (var onePointInNet2 in points2)
                    {
                        double subX = onePointInNet1.X - onePointInNet2.X;
                        double subY = onePointInNet1.Y - onePointInNet2.Y;
                        //是否可以用Compare函数
                        //onePointInNet1.Compare(onePointInNet2);
                        if ((subX * subX + subY * subY) < deviation)
                        {
                            mayMatchPointsInNet2.Add(onePointInNet2);
                        }
                    }
                    if (0 == mayMatchPointsInNet2.Count) //当情况为 1 ： 0
                    {


                    }
                    else if (1 == mayMatchPointsInNet2.Count) // 当情况为1 ：1
                    {
                        result1.Add(id, new PointF((float)onePointInNet1.X, (float)onePointInNet1.Y));
                        result2.Add(id, new PointF((float)mayMatchPointsInNet2[0].X, (float)mayMatchPointsInNet2[0].Y));
                        id++;
                    }
                    else // 当情况为 1：n  未考虑到m：n的情况  
                    {
                        //如何获取当前结点相连的边（即拓扑关系）
                        //var a = pGeomNet1.EdgeElement[onePointInNet1];
                        //var b = pGeomNet1.GeometryForEdgeEID[a];
                        //IPolyline line = b as IPolyline;
                        //var f = line.FromPoint;
                        //var t = line.ToPoint;

                    }
                }
            }

            return result;
        }

        //获取矩阵
        static Matrix GetMatrixFromPoint(Dictionary<int, System.Drawing.PointF> points)
        {
            if (points == null)
            {
                return null;
            }
            int siteCount = points.Count;
            Voronoi voroObject = new Voronoi();
            List<DelaunayTriangle> allTriangle = new List<DelaunayTriangle>();//delaunay三角形集合
            List<Site> sites = new List<Site>();
            List<Edge> trianglesEdgeList = new List<Edge>();//Delaunay三角形网所有边
            List<Edge> voronoiEdgeList = new List<Edge>();//vironoi图所有边
            List<Edge> voronoiRayEdgeList = new List<Edge>();//voroni图外围射线边
            foreach (var onePoint in points)
            {
                Site site = new Site(onePoint.Value.X, onePoint.Value.Y, onePoint.Key);
                sites.Add(site);
            }

            //将超级三角形的三点添加到三角形网中
            Site A = new Site(2500, -50000, -1);
            Site B = new Site(-50000, 4000, -1);
            Site C = new Site(50000, 4000, -1);
            DelaunayTriangle dt = new DelaunayTriangle(A, B, C);
            allTriangle.Add(dt);

            //构造Delaunay三角形网
            voroObject.setDelaunayTriangle(allTriangle, sites);

            //返回Delaunay三角形网所有边
            trianglesEdgeList = voroObject.returnEdgesofTriangleList(allTriangle);
            double[,] M = new double[siteCount, siteCount];
            for (int m = 0; m < siteCount; m++)
            {
                for (int n = 0; n < siteCount; n++)
                {
                    M[m, n] = 0;
                }
            }
            for (int i = 0; i < trianglesEdgeList.Count; i++)
            {
                int col = 0;
                int row = 0;
                row = trianglesEdgeList[i].a.ID;
                col = trianglesEdgeList[i].b.ID;
                if (row != -1 && col != -1)
                {
                    M[row, col] = 1;
                }
            }
            Matrix Mat = new Matrix(M);
            return Mat;
        }

        //m1 和 m2应当为维度相同的方阵
        static Matrix Disparity(Matrix m1, Matrix m2)
        {
            if ((m1.RowCount != m1.ColumnCount) || (m2.ColumnCount != m2.RowCount))
            {
                System.Console.WriteLine("不是方阵");
                return null;
            }
            if ((m1.ColumnCount != m2.ColumnCount))
            {
                System.Console.WriteLine("两个方阵维数不同");
                return null;
            }
            return m1 - m2;
        }

        //获取 j outline
        static int argmax(Matrix dM)
        {
            int arg = -1;
            if (dM == null)
            {
                return arg;
            }
            if (dM.RowCount != dM.ColumnCount)
            {
                System.Console.WriteLine("不是方阵");
                return arg;
            }
            int count = dM.RowCount;
            double max = 0.0;
            for (int i = 0; i < count; i++)
            {
                double sum = 0.0;
                for (int k = 0; k < count; k++)
                {
                    sum += dM[i, k];    //求和
                }
                if (sum > max)  //可以做一些优化吗？
                {
                    max = sum;  //和最大的行的行数 即为结果
                    arg = i;
                }
            }
            return arg;
        }

        static void Main(string[] args)
        {
            //arcGIS engine权限注册
            ESRI.ArcGIS.RuntimeManager.BindLicense(ProductCode.EngineOrDesktop);
            //test();
            //GetPoint("E:\\data\\test_1.gdb", "test_Net1");
            IGeometricNetwork pGNet1 = ReadNet("E:\\data\\arcgistest.gdb", "test", "test_Net2");
            IGeometricNetwork pGNet2 = ReadNet("E:\\data\\arcgistest.gdb", "test", "test_Net1");
            var a = OneToOneMatch(pGNet1, pGNet2);
            var pointsf1 = a[0];
            var pointsf2 = a[1];

            //List<Esri.FileGDB.Point> points1 = ReadFileToGetPoint("E:\\data\\test_1.gdb", "test\\test1");
            //List<Esri.FileGDB.Point> points2 = ReadFileToGetPoint("E:\\data\\test_1.gdb", "test\\test_Net1_Junctions");
            //Dictionary<int, System.Drawing.PointF> pointsf1 = Transform(points1);
            //Dictionary<int, System.Drawing.PointF> pointsf2 = Transform(points2);

            int j = -2;
            while (j != -1)
            {
                Matrix M1 = GetMatrixFromPoint(pointsf1);
                Matrix M2 = GetMatrixFromPoint(pointsf2);
                int count = pointsf1.Count;
                Matrix M = Disparity(M1, M2); //获取差异矩阵
                j = argmax(M);         //根据差异矩阵 获取 j outline
                if (j >= 0 && j < count)
                {
                    //pointsf1.Remove(j); //remove后记得将后面所有元素的key值前移
                    //pointsf2.Remove(j);
                    Removej(ref pointsf1, j);
                    Removej(ref pointsf2, j);
                }
                else
                {
                    break;
                }
            }

            //TestCreateFile(pointsf1);
            return;
        }


        //将FileGDB中的point转换 并为其编号
        //static Dictionary<int, System.Drawing.PointF> Transform(List<Esri.FileGDB.Point> points)
        //{
        //    Dictionary<int, System.Drawing.PointF> result = new Dictionary<int, PointF>();
        //    int id = 0;
        //    foreach (var one in points)
        //    {
        //        System.Drawing.PointF point = new System.Drawing.PointF(Convert.ToSingle(one.x), Convert.ToSingle(one.y));
        //        result.Add(id, point);
        //        id++;
        //    }
        //    return result;
        //}

        //remove 同时 将后面的点编号-1
        static void Removej(ref Dictionary<int, System.Drawing.PointF> points, int j)
        {
            Dictionary<int, System.Drawing.PointF> result = new Dictionary<int, PointF>();
            if (j >= 0 && j < points.Count)
            {
                points.Remove(j);
                foreach (var onePoint in points)
                {
                    if (onePoint.Key < j)
                    {
                        result.Add(onePoint.Key, onePoint.Value);
                    }
                    else
                    {
                        result.Add(onePoint.Key - 1, onePoint.Value);
                    }
                }
            }
            points = result;
        }

        //测试 将结果点集 生成 .shp文件
        //static void TestCreateFile(Dictionary<int, System.Drawing.PointF> points)
        //{
        //    Geodatabase geodatabase = Geodatabase.Open("E:\\data\\test_1.gdb");
        //    Esri.FileGDB.Table table = geodatabase.OpenTable("\\result");
        //    foreach (var one in points)
        //    {
        //        Esri.FileGDB.Row row = table.CreateRowObject();
        //        PointShapeBuffer geom = new PointShapeBuffer();
        //        geom.Setup(ShapeType.Point);
        //        Esri.FileGDB.Point point = new Esri.FileGDB.Point(one.Value.X, one.Value.Y);
        //        geom.point = point;
        //        row.SetGeometry(geom);
        //        table.Insert(row);
        //    }
        //    table.Close();
        //    geodatabase.Close();
        //    return;
        //}
    }
}

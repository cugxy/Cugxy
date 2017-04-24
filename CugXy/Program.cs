using System;
using System.Collections.Generic;
using System.IO;
using System.Drawing;
using VoronoiCSharpLib;
using LinearAlgebra;
using Esri.FileGDB;
using OSGeo.GDAL;
using OSGeo.OGR;
using ESRI.ArcGIS.Geodatabase;
using ESRI.ArcGIS.Controls;
using ESRI.ArcGIS.Geometry;
using ESRI.ArcGIS.DataSourcesGDB;
using ESRI.ArcGIS.NetworkAnalysis;
using ESRI.ArcGIS.Carto;
using ESRI.ArcGIS.Display;
using ESRI.ArcGIS.esriSystem;
using ESRI.ArcGIS;
using LinearAlgebra.VectorAlgebra;
using CugXy;

namespace CSharpTest
{

    class Program
    {
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
        }

        static List<Dictionary<int, PointF>> OneToOneMatch(IGeometricNetwork pGeomNet1, IGeometricNetwork pGeomNet2, double deviation = 3000.0)
        {
            List<Dictionary<int, PointF>> result = new List<Dictionary<int, PointF>>();
            Dictionary<int, PointF> result1 = new Dictionary<int, PointF>();
            Dictionary<int, PointF> result2 = new Dictionary<int, PointF>();
            if (pGeomNet1 == null)
            {
                Console.WriteLine("IGeometricNetwork 1 error");
                return null;
            }
            if (pGeomNet2 == null)
            {
                Console.WriteLine("IGeometricNetwork 2 error");
                return null;
            }
            int junctionCount1 = pGeomNet1.Network.JunctionCount;
            int junctionCount2 = pGeomNet2.Network.JunctionCount;
            DataStruct ds1;
            DataStruct ds2;
            if (junctionCount1 < junctionCount2)
            {
                ds1 = new DataStruct(pGeomNet1);
                ds2 = new DataStruct(pGeomNet2);
            }
            else
            {
                ds1 = new DataStruct(pGeomNet2);
                ds2 = new DataStruct(pGeomNet1);
            }
            int id = 0;
            foreach (var onePointInNet1 in ds1.MyPoints)
            {
                //pGeomNet2.JunctionElement[onePoint];
                //如何获取该点位置在第二张网络中附近的点呢
                //如果在第二张网络中遍历所有的点 判断距离 会不会太慢
                List<MyPoint> mayMatchPointsInNet2 = new List<MyPoint>();
                foreach (var onePointInNet2 in ds2.MyPoints)
                {
                    double subX = onePointInNet1.point.X - onePointInNet2.point.X;
                    double subY = onePointInNet1.point.Y - onePointInNet2.point.Y;
                    //是否可以用Compare函数
                    //onePointInNet1.Compare(onePointInNet2);
                    if ((subX * subX + subY * subY) < deviation)
                    {
                        mayMatchPointsInNet2.Add(onePointInNet2);
                    }
                }
                if (0 == mayMatchPointsInNet2.Count) //当情况为 1 ： 0
                {
                    //考虑是否欠妥
                    result1.Add(id, new PointF((float)onePointInNet1.point.X, (float)onePointInNet1.point.Y));
                    result2.Add(id, new PointF((float)mayMatchPointsInNet2[0].point.X, (float)mayMatchPointsInNet2[0].point.Y));
                    id++;
                }
                else if (1 == mayMatchPointsInNet2.Count) // 当情况为1 ：1
                {
                    result1.Add(id, new PointF((float)onePointInNet1.point.X, (float)onePointInNet1.point.Y));
                    result2.Add(id, new PointF((float)mayMatchPointsInNet2[0].point.X, (float)mayMatchPointsInNet2[0].point.Y));
                    id++;
                }
                else // 当情况为 1：n  未考虑到m：n的情况  
                {
                    //对于每一个待匹配结点
                    double maxSingn = 0.0;
                    int matchPoint = 0;
                    //foreach(var oneMatchPoint in mayMatchPointsInNet2)
                    for (int matchPointi = 0; matchPointi < mayMatchPointsInNet2.Count; matchPointi++)
                    {
                        var oneMatchPoint = mayMatchPointsInNet2[matchPointi];
                        //如何衡量 向量
                        //对于源结点中的每一条关联边
                        int num1 = onePointInNet1.lines.Count;//源结点的度
                        int num2 = oneMatchPoint.lines.Count;//待匹配结点的度
                        double[,] matchM = new double[num1, num2];
                        double Singn = 0.0;
                        int i = 0;
                        foreach (var oneEdge1 in onePointInNet1.lines)
                        {
                            int j = 0;
                            //构建向量
                            IPoint fPoint1 = oneEdge1.FromPoint;
                            IPoint tPoint1 = oneEdge1.ToPoint;
                            double length1 = oneEdge1.Length;
                            Vector XX = null;
                            if (onePointInNet1.point.X == fPoint1.X && onePointInNet1.point.Y == fPoint1.Y)
                            {
                                XX = new Vector(new double[2] { tPoint1.X - fPoint1.X, tPoint1.Y - fPoint1.Y }, VectorType.Row);
                            }
                            else
                            {
                                XX = new Vector(new double[2] { fPoint1.X - tPoint1.X, fPoint1.Y - tPoint1.Y }, VectorType.Row);
                            }
                            //对于待匹配结点的每一条关联边
                            foreach (var oneEdge2 in oneMatchPoint.lines)
                            {
                                IPoint fPoint2 = oneEdge2.FromPoint;
                                IPoint tPoint2 = oneEdge2.ToPoint;
                                Vector YY = null;
                                if (oneMatchPoint.point.X == fPoint2.X && oneMatchPoint.point.Y == fPoint2.Y)
                                {
                                    ICurve pLine;
                                    oneEdge2.GetSubcurve(0, length1, false, out pLine);
                                    IPoint fnPoint = pLine.FromPoint;
                                    IPoint tnPoint = pLine.ToPoint;
                                    YY = new Vector(new double[2] { tnPoint.X - fnPoint.X, tnPoint.Y - fnPoint.Y }, VectorType.Row);
                                    //按长度在边上取点后的向量
                                }
                                else
                                {
                                    ICurve pLine;
                                    double length2 = oneEdge2.Length;
                                    oneEdge2.GetSubcurve(length2, length2 - length1, false, out pLine);
                                    IPoint fnPoint = pLine.FromPoint;
                                    IPoint tnPoint = pLine.ToPoint;
                                    YY = new Vector(new double[2] { tnPoint.X - fnPoint.X, tnPoint.Y - fnPoint.Y }, VectorType.Row);
                                }//else 
                                 //计算向量匹配值 并填入矩阵中 
                                double matchValue = 0.0;
                                double XXLength = Math.Sqrt(XX.Elements[0] * XX.Elements[0] + XX.Elements[1] * XX.Elements[1]);
                                double YYLength = Math.Sqrt(YY.Elements[0] * YY.Elements[0] + YY.Elements[1] * YY.Elements[1]);
                                if (XXLength > YYLength)
                                {
                                    matchValue = (XX.Elements[0] * YY.Elements[0] + XX.Elements[1] * YY.Elements[1]) / (XXLength * XXLength);
                                }
                                else
                                {
                                    matchValue = (XX.Elements[0] * YY.Elements[0] + XX.Elements[1] * YY.Elements[1]) / (YYLength * YYLength);
                                }
                                matchM[i, j] = matchValue;
                                j++;
                            }//foreach  oneMatchPoint.Line
                            i++;
                        }//foreach  onePointInNet1.lines
                         //如果matchM[h, k]既是第 h 行最大值，又是第 k 列最大值，那么 h，k 之间存在匹配关系

                        for (int num1i = 0; num1i < num1; num1i++)
                        {
                            double tempMax = 0.0;
                            int tempRow = 0;
                            int tempCol = 0;
                            for (int h = 0; h < num1; h++)
                            {
                                for (int k = 0; k < num2; k++)
                                {
                                    if (matchM[h, k] > tempMax)
                                    {
                                        tempMax = matchM[h, k];
                                        tempRow = h;
                                        tempCol = k;
                                    }
                                }
                            }
                            //tempMax 一定大与等于 0 
                            Singn += tempMax;
                            for (int h = 0; h < num1; h++)
                            {
                                matchM[h, tempCol] = 0.0;
                            }
                            for (int k = 0; k < num2; k++)
                            {
                                matchM[tempRow, k] = 0.0;
                            }
                        }
                        if (Singn > maxSingn)
                        {
                            maxSingn = Singn;
                            matchPoint = matchPointi;
                        }
                    }// foreach  mayMatchPointsInNet2
                    result1.Add(id, new PointF((float)onePointInNet1.point.X, (float)onePointInNet1.point.Y));
                    result2.Add(id, new PointF((float)mayMatchPointsInNet2[matchPoint].point.X, (float)mayMatchPointsInNet2[matchPoint].point.Y));
                    id++;
                }//else 当为 1：n
            }//foreach MyPoints

            result.Add(result1);
            result.Add(result2);
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
            try
            {
                //arcGIS engine权限注册
                ESRI.ArcGIS.RuntimeManager.BindLicense(ProductCode.EngineOrDesktop);
                //test();
                //GetPoint("E:\\data\\test_1.gdb", "test_Net1");
                //string databasePath = string.Empty;
                //string dataSetName = string.Empty;
                //string net1Name = string.Empty;
                //string net2Name = string.Empty;
                //Console.WriteLine("请输入文件地理信息数据库路径，数据集名称，几何网络名称1，几何网络名称2");
                //databasePath = Console.ReadLine();
                //dataSetName = Console.ReadLine();
                //net1Name = Console.ReadLine();
                //net2Name = Console.ReadLine();
                //if (databasePath == string.Empty || dataSetName == string.Empty || net1Name == string.Empty || net2Name == string.Empty)
                //{
                //    Console.WriteLine("你的输入有误。");
                //    return;
                //}
                //IGeometricNetwork pGNet1 = ReadNet(databasePath, dataSetName, net1Name);
                //IGeometricNetwork pGNet2 = ReadNet(databasePath, dataSetName, net2Name);
                IGeometricNetwork pGNet1 = ReadNet("E:\\data\\arcgistest.gdb", "test", "test_Net2");
                IGeometricNetwork pGNet2 = ReadNet("E:\\data\\arcgistest.gdb", "test", "test_Net1");

                List<Dictionary<int, PointF>> pointsList = OneToOneMatch(pGNet1, pGNet2);
                if (pointsList == null)
                {
                    Console.WriteLine("One to One Match error. without result");
                    return;
                }

                //List<Esri.FileGDB.Point> points1 = ReadFileToGetPoint("E:\\data\\test_1.gdb", "test\\test1");
                //List<Esri.FileGDB.Point> points2 = ReadFileToGetPoint("E:\\data\\test_1.gdb", "test\\test_Net1_Junctions");
                //Dictionary<int, System.Drawing.PointF> pointsf1 = Transform(points1);
                //Dictionary<int, System.Drawing.PointF> pointsf2 = Transform(points2);

                Dictionary<int, System.Drawing.PointF> pointsf1 = pointsList[0];
                Dictionary<int, System.Drawing.PointF> pointsf2 = pointsList[1];
                if (pointsf1 == null || pointsf2 == null)
                {
                    Console.WriteLine("One to One Match error.");
                    return;
                }

                TestCreateFile(pointsf1, "E:\\data\\Result.gdb", "firstResult1");
                TestCreateFile(pointsf2, "E:\\data\\Result.gdb", "firstResult2");

                int j = -2;
                while (j != -1)
                {
                    Matrix M1 = GetMatrixFromPoint(pointsf1);
                    Matrix M2 = GetMatrixFromPoint(pointsf2);
                    int count = pointsf1.Count;
                    Matrix M = Disparity(M1, M2); //获取差异矩阵
                    if (M == null)
                    {
                        Console.WriteLine("fiald to get M.");
                        return;
                    }
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

                TestCreateFile(pointsf1, "E:\\data\\Result.gdb", "result");
                return;
            }
            catch (Exception e)
            {
                Console.WriteLine("error" + e.Message);
                return;
            }
        }


        //将FileGDB中的point转换 并为其编号
        static Dictionary<int, System.Drawing.PointF> Transform(List<Esri.FileGDB.Point> points)
        {
            Dictionary<int, System.Drawing.PointF> result = new Dictionary<int, PointF>();
            int id = 0;
            foreach (var one in points)
            {
                System.Drawing.PointF point = new System.Drawing.PointF(Convert.ToSingle(one.x), Convert.ToSingle(one.y));
                result.Add(id, point);
                id++;
            }
            return result;
        }

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
        static void TestCreateFile(Dictionary<int, System.Drawing.PointF> points, string filePath, string tableName)
        {
            Geodatabase geodatabase = Geodatabase.Open(filePath);
            Esri.FileGDB.Table table = geodatabase.OpenTable(tableName);
            foreach (var one in points)
            {
                Esri.FileGDB.Row row = table.CreateRowObject();
                PointShapeBuffer geom = new PointShapeBuffer();
                geom.Setup(ShapeType.Point);
                Esri.FileGDB.Point point = new Esri.FileGDB.Point(one.Value.X, one.Value.Y);
                geom.point = point;
                row.SetGeometry(geom);
                table.Insert(row);
            }
            table.Close();
            geodatabase.Close();
            return;
        }
    }
}

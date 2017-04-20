
using System;
using System.Collections.Generic;
using System.Drawing;
namespace VoronoiCSharpLib
{

    //点
    public class Site
	{
		public double x, y;
        public int ID;
        public Site()
        { }
        public Site(double x, double y, int ID)
		{
            this.x = x;
            this.y = y;
            this.ID = ID;
		}
	}


    //三角形的边
    public class Edge
	{
        public Site a, b;
        public Edge(Site a, Site b)
        {
            this.a = a;
            this.b = b;
        }
	}

    //自定义排序规则
    public class SiteSorterXY : IComparer<Site>
	{
        public int Compare(Site p1, Site p2)
		{
			if ( p1.x > p2.x ) return 1;
            if (p1.x < p2.x) return -1;
			return 0;
		}
	}

    public class DelaunayTriangle
    {
        Voronoi voronoi = new Voronoi();
        public Site site1, site2, site3;//三角形三点
        public Site centerPoint;//外界圆圆心
        public double radius;//外接圆半径
        public List<DelaunayTriangle> adjoinTriangle;//邻接三角形 

        public DelaunayTriangle(Site site1,Site site2,Site site3)
        {
            centerPoint = new Site();
            this.site1 = site1;
            this.site2 = site2;
            this.site3 = site3;
            //构造外接圆圆心以及半径
            voronoi.circle_center(centerPoint, site1, site2,site3,ref radius);
        }
    }


}

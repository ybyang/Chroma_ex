#include "inline_ex/inline_Ob.h"
#include "chroma.h"
#include <math.h>
#include <complex>
#include <iostream>
#include <string>


namespace Chroma 
{
    
    namespace InlineObEnv 
    {
	//Name of the measurement to be called in the XML input file
	const std::string name = "GMF_O_b";
	
	//This function is used with the factory thing
	AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
						const std::string& path) 
	{
	    //Create new instance of measurement class using params
	    //from the passed xml file
	    return new InlineMyMeas(InlineObParams(xml_in, path));
	}
		
	// Local registration flag
	namespace {
	    bool registered = false;
	}
	
	// Function to register all the factories
	bool registerAll() 
	{
	    bool success = true; 
	    if (! registered)
	    {
		success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
		QDPIO::cout << "Registering " << name << " " << success << std::endl;
		registered = true;
	    }
	    return success;
	}
	
	
	
	/*** Implementation of Parameter functions ***/
	
	//Set default parameters
	InlineObParams::InlineObParams() { frequency = 0; radius = 0; srcs.resize(1); }
	
	//Read parameters in from xml file
	InlineObParams::InlineObParams(XMLReader& xml_in, const std::string& path) 
	{
	    try 
	    {
		XMLReader paramtop(xml_in, path);
		
		if (paramtop.count("Frequency") == 1)
		    read(paramtop, "Frequency", frequency);
		else
		    frequency = 1;
		
		read(paramtop, "NamedObject", named_obj);

		//Read in the starting position time (ie the time
		//of the first source) and the time between each
		//source. This assumes that each source is equally
		//spaced in time. See branch First-O_b for code
		//that drops this assumtion and takes sources
		//as the arguments (loc & start/stop times)
		read(paramtop, "Multi_Src", srcs);

		if(paramtop.count("radius") == 1)
		    read(paramtop, "radius", radius);
		else
		    radius = 0;

		
		// Possible alternate XML file pattern
		if (paramtop.count("xml_file") != 0) 
		{
		    read(paramtop, "xml_file", xml_file);
		} else
		    xml_file = "";
		
	    }
	    catch(const std::string& e) 
	    {
		QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
		QDP_abort(1);
	    }
	}
	
	
	// Write loaded params to xml file
	void InlineObParams::write(XMLWriter& xml_out, const std::string& path) 
	{
	    push(xml_out, path);
	    
	    // write our all params
	    QDP::write(xml_out, "Multi_Src", srcs);
	    //QDP::write(xml_out, "Named_Object", named_obj);
	    QDP::write(xml_out, "radius", radius);

		  
	    
	    if(xml_file != "")
		QDP::write(xml_out, "xml_file", xml_file);
	    pop(xml_out);
	}
    } //End namespace InlineObEnv

    /*** Inline Measurement function implimentation ***/
    // Function call
    void InlineMyMeas::operator()(unsigned long update_no,
				  XMLWriter& xml_out) 
    {
	// If xml file not empty, then use alternate
	if (params.xml_file != "")
	{
	    std::string xml_file = makeXMLFileName(params.xml_file, update_no);

	    push(xml_out, "GMF_O_b");
	    write(xml_out, "update_no", update_no);
	    write(xml_out, "xml_file", xml_file);
	    pop(xml_out);
	    
	    XMLFileWriter xml(xml_file);
	    func(update_no, xml);
	}
	else
	{
	    func(update_no, xml_out);
	}
    }
    
    
    /*** Measurement code stars here ***/
    void InlineMyMeas::func(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
	START_CODE();
	
	QDPIO::cout << InlineObEnv::name << ": Begining" << std::endl;
	
	StopWatch snoop;
	snoop.reset();
	snoop.start();

	//Print out boilerplate stuff to xml
	push(xml_out, "GMF_O_b");
	write(xml_out, "update_no", update_no);

	//Write out the input
	params.write(xml_out, "Input");


	/** Calculate the two dimensional plaquettes **/
	
	//Get link matrix
	multi1d<LatticeColorMatrix> u;
	u = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

	//Create variables to store the plaquette and the average
	//trace of the plaquette
	multi3d<LatticeColorMatrix> plane_plaq;
	multi2d<Double> tr_plane_plaq;
	plane_plaq.resize(4,Nd,Nd); //Quadrent, plane
	tr_plane_plaq.resize(Nd,Nd);
	Double w_plaq;
	Double s_plaq;
	Double t_plaq;

	/* Calculate the plaquette and the average trace of
	 * the plaquette
	 */
	for(int mu = 0; mu < Nd; mu++)
	{
	    for(int nu = mu+1; nu < Nd; nu++)
	    {
		LatticeColorMatrix tmp, tmp2, tmp3;
		//Do first quadrent
		tmp = shift(u[nu], FORWARD, mu);
		tmp2 = u[mu] * tmp;
		tmp = shift(u[mu], FORWARD, nu);
		tmp3 = u[nu] * tmp;

		//Record the plaquette in a cross section
		plane_plaq[0][nu][mu] = tmp2*adj(tmp3);

		//Do second quadrent
		tmp = shift(shift(u[mu],FORWARD, nu), BACKWARD, mu);
		tmp2 = u[nu] * adj(tmp);
		tmp = shift(u[nu], BACKWARD, mu);
		tmp3 = adj(tmp) * shift(u[mu], BACKWARD, nu);

		plane_plaq[1][nu][mu] = tmp2*tmp3;

		//Do third quad
		tmp = shift(u[mu], BACKWARD, mu);
		tmp2 = shift(shift(u[nu], BACKWARD, nu), BACKWARD, mu);
		tmp2 = adj(tmp) * adj(tmp2);
		tmp = shift(shift(u[mu], BACKWARD, nu), BACKWARD, mu);
		tmp3 = tmp * shift(u[nu], BACKWARD, nu);

		plane_plaq[2][nu][mu] = tmp2*tmp3;

		//Do fourth quad
		tmp = shift(u[nu], BACKWARD, nu);
		tmp2 = adj(tmp) * shift(u[mu], BACKWARD, nu);
		tmp3 = shift(shift(u[nu], FORWARD, mu), BACKWARD, nu) * adj(u[mu]);

		plane_plaq[3][nu][mu] = tmp2*tmp3;
		
		tr_plane_plaq[nu][mu] = sum(real(trace(plane_plaq[0][nu][mu])));
		
		//Normalize the plane
		tr_plane_plaq[nu][mu] /= Double(Layout::vol() * Nc);
		plane_plaq[0][mu][nu] = plane_plaq[0][nu][mu]; //symmetric
		tr_plane_plaq[mu][nu] = tr_plane_plaq[nu][mu]; //symmetric
		
		w_plaq += tr_plane_plaq[nu][mu];

		//record the time/space plaq
		if(nu == tDir())
		    t_plaq += tr_plane_plaq[nu][mu];
		else
		    s_plaq += tr_plane_plaq[nu][mu];

	    } //end nu loop
	} //Found normalized plane and unorm plaq: end mu loop

	//Normalize the average plaquette for space/time/whole
	w_plaq *= 2 / Double(Nd*(Nd-1));
	t_plaq /= (Nd - 1);
	if(Nd > 2)
	    s_plaq *= 2 / Double((Nd-2)*(Nd-1));
	

	/** Record the results in the xml file **/
	write(xml_out, "w_plaq", w_plaq);
	write(xml_out, "s_plaq", s_plaq);
	write(xml_out, "t_plaq", t_plaq);

	/** Write plane plaq to xml file **/
	
	for(int mu = 0; mu < Nd; mu++)
	{
	    for(int nu = mu+1; nu < Nd; nu++)
	    {
		write(xml_out, "plane_plaq_" + std::to_string(mu) +
		      std::to_string(nu), tr_plane_plaq[mu][nu]);
	    }
	}
	
	QDPIO::cout << "Finding F" << std::endl;
	/** Find F_{n,munu} **/
	multi2d<LatticeColorMatrix> F;
	F.resize(Nd,Nd);
	for(int mu = 0; mu < Nd; mu++)
	{
	    for(int nu = mu+1; nu < Nd; nu++)
	    {
		for(int i = 0; i < 4; i++)
		    F[mu][nu] += plane_plaq[i][mu][nu] - adj(plane_plaq[i][mu][nu]);

		//TODO: Add factor of 1/(8 i g a^2)
		F[mu][nu] *= 1/8.0;
		F[nu][mu] = F[mu][nu]; //Symmetric
	    }
	}
	QDPIO::cout << "Finding E/B" << std::endl;
	/** Calculate E and B**/
	multi1d<LatticeColorMatrix> E,B;
	E.resize(3);
	B.resize(3);
	for(int i = 0; i < 3; i++)
	  {
	    E[i] = F[3][i];
	    for(int j = 0; j < 3; j++)
	      for(int k = 0; k<i; k++)
		B[i] += leviCh(i,j,k)*F[j][k];
	  }


	std::vector<Double> VecOa;
	QDPIO::cout << "Finding Oa" << std::endl;
	LatticeColorMatrix Oa;
	for(int i = 0; i < 3;i++)
	  Oa += E[i]*E[i]-B[i]*B[i];
		      
	QDPIO::cout << "Finding Tr(Oa)" << std::endl;
	for(int t = 0; t < Layout::lattSize()[3]; t++)
	  {
	    Double a = 0;
	    multi1d<int> tCoords;
	    tCoords.resize(Nd);
	    tCoords[3] = t;
	    for(int x = 0; x < Layout::lattSize()[0];x++)
	      for(int y = 0; y < Layout::lattSize()[1];y++)
		for(int z = 0; z < Layout::lattSize()[2];z++)
		{
		  tCoords[0] = x; tCoords[1] = y; tCoords[2] = z;
		  a += real(trace(peekSite(Oa, tCoords)));
		}

	    VecOa.push_back(a);
	  }

	write(xml_out, "Oa", VecOa);
       
	
	/*** Calculating O_b ***/
	std::vector<Double> O_b;
	std::vector<Double> E2;
	std::vector<Double> B2;

        for(int i = 0; i < params.srcs.size(); i++)
	{
	    //Start at t_start
	    QDPIO::cout << "Processing src " << i << " from "
			<< params.srcs[i].t_start << " to "
			<< params.srcs[i].t_end << std::endl;
	    getO_b(O_b, E2, B2, params.srcs[i], params.radius, plane_plaq[0]);
	}
	    
	//Write O_b out to the xml file
	write(xml_out, "O_b", O_b);
	write(xml_out, "E2", E2);
	write(xml_out, "B2", B2);

	pop(xml_out);
		
	snoop.stop();
	QDPIO::cout << InlineObEnv::name << ": total time = "
		    << snoop.getTimeInSeconds() 
		    << " secs" << std::endl;
	
	QDPIO::cout << InlineObEnv::name << ": ran successfully" << std::endl;
	
	END_CODE();
    } //End func()

    ColorMatrix InlineMyMeas::get_G(const multi1d<int>& coords, int mu, int nu,
				    const multi2d<LatticeColorMatrix>& P)
    {
	//TODO impliment boundry conditions
	ColorMatrix G;
	//Should pull in the four plaquettes
	for(int i = 0; i <= 1; i++)
	    for(int j = 0; j <= 1; j++)
	    {
		multi1d<int> tCoords = coords;
		tCoords[mu] -= i;
		tCoords[nu] -= j;
		G += peekSite(P[mu][nu], tCoords);
		QDPIO::cout << "Getting P_{" << mu << nu << "} at ( ";
		for(int x = 0; x < Nd; x++)
		    QDPIO::cout << tCoords[x] << " ";
		QDPIO::cout << ")" << std::endl;
	    }
	
	return G;
    }


    //Code to sum over all spaceial positions of O_b and return a vector O_b(t) = O_b[t]. 
    void InlineMyMeas::getO_b(std::vector<Double>& vecOb,
			      std::vector<Double>& vecE,
			      std::vector<Double>& vecB,
			      const  InlineObEnv::InlineObParams::Src_t src,
			      const int radius,
			      const multi2d<LatticeColorMatrix>& plane_plaq)
    {
	//Constants for finding O_b
	multi1d<int> t_coords; //Coords to find O_b at
	t_coords.resize(Nd);
	Double Beta = 1;
	Double a = 1;
	
	for(int t = src.t_start; t != src.t_end+1; t = (t+1)% Layout::lattSize()[3] )
	{
	    //QDPIO::cout << "Processing t=" << t << std::endl;
	    Double O_b = 0;
	    Double E = 0;
	    Double B = 0;
	    t_coords[3] = t;

	    /** Sum over all space **/
	    for(int x = 0; x < Layout::lattSize()[0]; x++)
	    {
		t_coords[0] = x;
		for(int y = 0; y < Layout::lattSize()[1]; y++)
		{
		    t_coords[1] = y;
		    for(int z = 0; z < Layout::lattSize()[2]; z++)
		    {
			t_coords[2] = z;
			if(validLocation(t_coords, src.srcLoc, radius))
			{
			    Double e,b;
			    O_b += getO_b(t_coords, plane_plaq, e, b);
			    E += e;
			    B += b;
			}
		    } 
		}
	    } //End sum over space

	    /* add scale factors to O_b and push back to
	       resultant vector 
	    */
	    //TODO: Should this be scaled by 1/Nc like the normalized plaquettes are?
	    //O_b *= -1 * Double(4)/9.0 * Beta/a * Double(1)/Nc;
	    vecOb.push_back(O_b);
	    vecE.push_back(E);
	    vecB.push_back(B);
	} //end loop through time
    } //end getO_b
    
    //Code to calculate O_b at coords
    Double InlineMyMeas::getO_b(const multi1d<int>& t_coords,
				const multi2d<LatticeColorMatrix>& plane_plaq,
				Double& E, Double& B)
    {
	E = 0;
	B = 0;
	
	for(int mu = 0; mu < Nd; mu++)
	{
	    if(mu != tDir())
	    {
		//Do first half of sum
		E += real(trace(
				peekSite(plane_plaq[mu][tDir()],
					 t_coords)));
		//do second half of sum
		for(int nu = 0; nu < mu; nu++)
		{
		    if(nu != tDir())
		    {
			B += real(trace(
					peekSite(plane_plaq[nu][mu],
						 t_coords)));
		    }
		}//end second half
	    } //end check for tDir
	} //end 1st half
        return E-B;
    } //end getO_b

    bool InlineMyMeas::validLocation(const multi1d<int>& t_coords,
				     const multi1d<int>& t_src,
				     int R)
    {
	if (R == 0) return true;
	
	int dist = 0;
	for(int i = 0; i < 3; i++)
	{
	    int dx = std::abs(t_coords[0] - t_src[0]);
	    int dimSize = Layout::lattSize()[i];
	    dx = (dx > dimSize - dx) ? dimSize - dx : dx;
	    dist += dx*dx;
	}
	
	return dist <= R*R;
	    
    }
};

#pragma once

namespace CFDMesh {
namespace Equation {

template<typename Config>
struct Equation {
	using real = typename Config::real;
	
	virtual ~Equation() {}
	
	//TODO abstract away the arguments ...
	//or just turn all of Equation into a template parameter, and make everything to do with Equation compile-time
	//virtual std::pair<real, real> calcLambdaMinMax(vec normal, Prim W, real Cs) = 0;
};

}
}

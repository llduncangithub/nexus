

THREE.NXZLoader = function(manager) {
    this.manager = (manager !== undefined) ? manager : THREE.DefaultLoadingManager;
};


THREE.NXZLoader.prototype = {

    constructor: THREE.PLYLoader,

    load: function(url, onLoad, onProgress, onError) {
        const scope = this;
        const loader = new THREE.FileLoader(scope.manager);
        loader.setPath(this.path);
        loader.setResponseType('arraybuffer');
        loader.load(url, function(blob) {
            onLoad(scope.decodeNxz(blob));
        }, onProgress, onError);
    },

    decodeNxz: function(buffer) {

		var now = performance.now();
		for(var i = 0; i < 10; i++) {
	        var decoder = new NxzDecoder(buffer);
	        var mesh = decoder.decode();	
		}
		console.log("MT/s:", (10*mesh.nface/1000000)/((performance.now() - now)/1000));

//Mesh is an an array made like this.
/*        mesh = {
            index:    new Uint32Array(nface*3),
            position: new Float32Array(nvert*3),
            normal:   new Float32Array(nvert*3),
            uv:       new Float32Array(nvert*2),
            color:    new Float32Array(nvert)
        };  */

        var geometry = new THREE.BufferGeometry();
        if (mesh.nface) 
          geometry.setIndex(new THREE.BufferAttribute(mesh.index, 1));

		geometry.addAttribute('position', new THREE.BufferAttribute(mesh.position, 3));
		if(mesh.color)
			geometry.addAttribute('color', new THREE.BufferAttribute(mesh.color, 3, true));
        if (mesh.normal)
          geometry.addAttribute('normal', new THREE.BufferAttribute(mesh.normal, 3, true));
        if (mesh.uv)
            geometry.addAttribute('uv', new THREE.BufferAttribute(mesh.uv, 2));
		return geometry;
    }
};

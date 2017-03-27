

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

        var decoder = new NxzDecoder(buffer);
        var mesh = decoder.decode();

//Mesh is an an array made like this.
/*        mesh = {
            index:    new Uint32Array(nface*3),
            position: new Float32Array(nvert*3),
            normal:   new Float32Array(nvert*3),
            uv:       new Float32Array(nvert*2),
            color:    new Float32Array(nvert)
        };  */

        const geometry = new THREE.BufferGeometry();
        if (nface) {
          geometry.setIndex(new(mesh.position.length/3 > 65535 ?
                THREE.Uint32BufferAttribute : THREE.Uint16BufferAttribute)
              (mesh.index, 1));
        }
        geometry.addAttribute('position', new THREE.Float32BufferAttribute(mesh.position, 3));
		if(mesh.color)
	        geometry.addAttribute('color', new THREE.Float32BufferAttribute(mesh.color, 3));
        if (mesh.normal)
          geometry.addAttribute('normal', new THREE.Float32BufferAttribute(mesh.normals, 3));
        
        if (mesh.uv)
            geometry.addAttribute('uv', new THREE.Float32BufferAttribute(mesh.uvs, 2));
		return geometry;
    }
};

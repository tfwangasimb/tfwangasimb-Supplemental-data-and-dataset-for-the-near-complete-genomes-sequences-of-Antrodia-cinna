#!/usr/bin/env node
// @ts-check
//
// input table, required headers: [GeneID, TranscriptID, Contig, Start, Stop, InterPro, GO Terms, antiSMASH, CAZyme]
// node caz.js --i=aaa.txt --min-caz=3 --max-ns=2

const fs = require("fs");


const min_caz = 3;
const max_ns = 2;

/**
 * transcription factors list
 */
const TFs_list = tsv_parse("TFs.txt").map(a => a[0]);

/**
 * isTC
 * @example
 * RawGene["InterPro"].indexOf("transporter") >= 0;
 * RawGene["GO Terms"].indexOf("transporter") >= 0;
 */
const Transporter = "transporter";
const input_annotation_table_CAZyme_columnName = "dbCAN";

class Gene {
	constructor() {
		this.GeneID = "";
		this.Contig = "";

		/** @type {number|string} */
		this.Start = 0;

		/** @type {number|string} */
		this.Stop = 0;
		
		/** @type {"CAZyme"|"TF"|"TC"|"NS"} */
		this.type = "NS";// Non-signature gene

		this.antiSMASH = null;

		/** @type {string} */
		this.smId = null;
	}
}

class RawGene extends Gene {
	GeneID = "";
	Contig = "";
	// Start = "";
	// Stop = "";
	Strand = "";
	InterPro = "";// [TFs, TC]
	GO_Terms = "";// [TC]
	antiSMASH = "";
	
	CAZyme =  "";
	
	TranscriptID = "";
	Feature = "";
}

const header_mapping = {
	"GeneID": "GeneID",
	"Contig": "Contig",
	"Start": "Start",
	"Stop": "Stop",
	"Strand": "Strand",
	
	"InterPro": "InterPro",// [TFs, TC]
	"GO Terms": "GO_Terms",// [TC]

	[input_annotation_table_CAZyme_columnName]: "CAZyme",

	"antiSMASH": "antiSMASH",
	
	"TranscriptID": "TranscriptID",
	"Feature": "Feature",
};

/**
 * @param {{ [head:string]: string; }} obj
 * @param {string} rawHeader
 * @param {string} val
 */
function table_transform(obj, rawHeader, val) {
	const h = header_mapping[rawHeader];
	if (h) {
		obj[h] = val;
	}
	return obj;
}

/**
 * @template T
 * @param {T[]} list
 * @param {string} key
 * @returns {T[][]}
 */
function groupBy(list, key) {
	const g = {
	};
	list.forEach(a => {
		if (!g[a[key]]) {
			g[a[key]] = [];
		}
		g[a[key]].push(a);
	});
	return Object.values(g);
}

/**
 * @param {RawGene[]} _rows input
 */
function mergeGeneIsoform(_rows) {
	/** @type {Map<string, RawGene[]>} */
	const map = new Map();

	console.log("rows.length", _rows.length);

	_rows.forEach(a => {
		if (a.Start > a.Stop) {
			[a.Start, a.Stop] = [a.Stop, a.Start];
		}
	});

	const rows = groupBy(_rows, "Contig").reduce((flat, g) => {
		const s = g.sort((a, b) => {
			const d1 = Number(a.Start) - Number(b.Start);
			const d2 = Number(a.Stop) - Number(b.Stop);
			return d1 != 0 ? d1 : d2;
		});
		return [].concat(flat, s);
	}, []);

	for (let i = 0; i < rows.length; ++i) {
		let cols = rows[i];
		const { GeneID, Contig } = cols;
		const k = [GeneID, Contig].join("_");
		const isoform_list = map.get(k) || [];
		isoform_list.push(cols);
		map.set(k, isoform_list);
	}

	console.log("map.size", map.size);

	/** @type {RawGene[]} */
	const new_list = [];

	map.forEach(function (isoform_list, k) {
		if (isoform_list.length > 1) {
			
			let start = Infinity;
			let end = -Infinity;
			isoform_list.forEach(cols => {
				let { Start, Stop } = cols;
				
				// sort
				[Start, Stop] = [Start, Stop].sort((a, b) => Number(a) - Number(b));

				start = Math.min(start, Number(Start));
				end = Math.max(end, Number(Stop));
			});

			// const { GeneID, Contig } = isoform_list[0];
			
			const merged = Object.assign(new RawGene(), isoform_list[0]);
			// merged.GeneID = merged.GeneID;
			// merged.Contig = merged.Contig;
			merged.Start = String(start);
			merged.Stop = String(end);
			//
			merged.InterPro = distinct_notEmpty(isoform_list.map(a => a.InterPro)).join(";");
			merged.GO_Terms = distinct_notEmpty(isoform_list.map(a => a.GO_Terms)).join(";");
			merged.CAZyme = distinct_notEmpty(isoform_list.map(a => a.CAZyme)).join(";");
			merged.TranscriptID = distinct_notEmpty(isoform_list.map(a => a.TranscriptID)).join(";");

			function distinct_notEmpty(list) {
				const s = new Set(list);
				return [...s].filter(a => a);
			}

			new_list.push(merged);
		}
		else if (isoform_list.length == 1) {
			new_list.push(isoform_list[0]);
		}
		else {
			throw new Error("err isoform: " + k);
		}
	});

	return new_list;
}

/**
 * @param {string} saveName
 * @param {string} fileName - funannotate annotations
 * @param {(obj: { [head:string]: string; }, rawHeader: string, val: any) => { [head:string]: string; }} table_transform
 * @param {number} min_caz
 * @param {number} max_ns
 */
function CAZyme_cluster(saveName, fileName, table_transform, min_caz, max_ns) {
	const gene_list = (function () {
		/** @type {any[]} */
		const _list = table_to_object_list(tsv_parse(fs.readFileSync(fileName).toString()), 0, { transform: table_transform, });
		/** @type {RawGene[]} */
		const list = _list;
		return mergeGeneIsoform(list);
	})();

	console.log(...(new Set(gene_list.map(a => a.Contig))));

	console.log(fileName, gene_list.length);

	const fstream = fs.createWriteStream(saveName);

	// SM cluster
	{
		/**
		 * @type {{ [smId: string]: Gene[]; }}
		 */
		const sm_geneList = {
		};

		const SM_cluster = gene_list.forEach(_gene => {
			/** @type {Gene} */
			const gene = (_gene);
			if (gene.antiSMASH) {
				const smId = gene.antiSMASH.match(/\d+/)[0];
				gene.smId = smId;
				sm_geneList[smId] = sm_geneList[smId] || [];
				sm_geneList[smId].push(gene);
			}
		});

		const smList = Object.values(sm_geneList);
		
		const sm_region_list = smList.map((gene_list, sm_idx) => {
			const start = Math.min(...gene_list.map(b => Number(b.Start)));
			const end = Math.max(...gene_list.map(b => Number(b.Stop)));
			
			const nChr = get_nChr(gene_list);
			return [
				"cluster",
				nChr,
				gene_list[0].smId, //`${nChr}.${smId}`,
				start,
				end,
				"#FF0000",
				gene_list.length,
				"#FF0000",
				gene_list.map(a => a.GeneID).join(";"),
			].join("\t");
		});

		if (sm_region_list.length) {
			fstream.write(sm_region_list.join("\n"));
			fstream.write("\n");
		}
	}

	// CAZyme cluster
	{
		const { caz_cluster_list, result_table } = _CAZyme_cluster(gene_list, min_caz, max_ns);

		/**
		 * @type {{ [nChr: number]: number; }}
		 */
		const chr_ca_last_SN = {
		};

		const caz_region_list = caz_cluster_list.map(gene_list => {
			const start = Math.min(...gene_list.map(b => Number(b.Start)));
			const end = Math.max(...gene_list.map(b => Number(b.Stop)));
			
			const nChr = get_nChr(gene_list);
			
			chr_ca_last_SN[nChr] = (chr_ca_last_SN[nChr] | 0) + 1;
			const caId = chr_ca_last_SN[nChr];//caIdx + 1;
			
			// type	chr	name	start	end	color	desc
			return [
				"CAZyme-cluster",
				nChr,
				`${nChr}.${caId}`,// `CA${nChr}.${caId}`,
				start,
				end,
				"#00EEFF",
				gene_list.length,
				"#0000FF",
				gene_list.map(a => a.GeneID).join(";"),
			].join("\t");
		});

		fstream.write(caz_region_list.join("\n"));
		fstream.write("\n");
	}

	// cluster gene
	{
		const caz_gene_list = gene_list.filter(gene => gene["CAZyme"]).map(gene => {
			const nChr = get_nChr([gene]);
			return [
				"gene",
				nChr,
				gene.CAZyme,
				gene.Start,
				gene.Stop,
				"#7F7F7F",
				"",
				"#000000",
			].join("\t");
		});
		
		fstream.write(caz_gene_list.join("\n"));
		fstream.write("\n");
	}
}

/**
 * @param {(Gene|RawGene)[]} gene_list
 */
function get_nChr(gene_list) {
	if (!gene_list.every(v => v.Contig == gene_list[0].Contig)) {
		console.error(...(new Set(gene_list.map(a => a.Contig))));
		throw new Error("if (gene_list.every(v => v.Contig == gene_list[0].Contig)) {");
	}

	const m = gene_list[0].Contig.match(/\d+$/);
	if (m) {
		return m[0];
	}
	else {
		return gene_list[0].Contig;
	}
}

/**
 *  @param {RawGene[]} gene_list
 *  @param {number} min_caz
 *  @param {number} max_ns
 */
function _CAZyme_cluster(gene_list, min_caz, max_ns) {
	/** @type {Gene[][]} */
	const caz_cluster_list = [];

	/** @type {Gene[]} */
	let current_cluster = [];

	const counter = {
		A: 0,
		F: 0,
		C: 0,
		N: 0,
	};

	const result_table = gene_list.map(cols => {
		const { GeneID, TranscriptID, Contig, Start, Stop, InterPro, GO_Terms, CAZyme } = cols;
		
		const last = current_cluster[current_cluster.length - 1];
		if (last != null && last.Contig != Contig) {
			end_cluster();
		}

		function push_gene(type) {
			const gene = new Gene();
			gene.GeneID = TranscriptID || GeneID;
			// gene.TranscriptID = TranscriptID;
			gene.Contig = Contig;
			gene.Start = Number(Start);
			gene.Stop = Number(Stop);
			gene.type = type;
			current_cluster.push(gene);
		}

		const isCAZyme = !!CAZyme;
		const isTF = InterPro != null && TFs_list.some(tf => InterPro.indexOf(tf) >= 0);
		const isTC = (
			(InterPro != null && InterPro.indexOf(Transporter) >= 0) ||
			(GO_Terms != null && GO_Terms.indexOf(Transporter) >= 0)
		);

		if (isCAZyme) {
			push_gene("CAZyme")
			counter["A"] += 1;
			counter["N"] = Math.max(0, counter["N"] - 1);
			// return "A";
			// return "CAZyme";
		}
		else if (isTF) {
			push_gene("TF");
			counter["F"] += 1;
			counter["N"] = Math.max(0, counter["N"] - 1);
			// return "F";
			// return "TF";
		}
		else if (isTC) {
			push_gene("TC");
			counter["C"] += 1;
			counter["N"] = Math.max(0, counter["N"] - 1);
			// return "C";
			// return "TC";
		}
		else {
			counter["N"] += 1;

			if (counter["N"] <= max_ns) {
				if (current_cluster.length) {
					push_gene("NS");
				}
			}
			else {
				end_cluster();
			}
			// return "N";
			// return "Non-signature";
		}

		return [
			GeneID,
			isCAZyme ? CAZyme : "",
			isTF ? "TF" : "",
			isTC ? "TC" : "",
			TranscriptID,
		];
	});
	return { caz_cluster_list, result_table };

	function end_cluster() {
		if (current_cluster.length > 0 &&
			current_cluster[current_cluster.length - 1].type == "NS"
		) {
			current_cluster.pop(); // remove last Non-signature gene
		}

		if (CAZyme_definition(counter, min_caz)) {
			caz_cluster_list.push(current_cluster);
		}
		current_cluster = []; // new
		counter["A"] = 0;
		counter["F"] = 0;
		counter["C"] = 0;
		counter["N"] = 0;
	}
}

/**
 * @param {{ A: number; F: number; C: number; N: number; }} counter
 * @param {number} min_caz
 * @returns {boolean}
 */
function CAZyme_definition(counter, min_caz) {
	const r = [
		(
			counter["A"] >= min_caz &&
			counter["F"] == 0 &&
			counter["C"] == 0
		),
		(
			counter["A"] >= (min_caz - 1) &&
			(counter["F"] + counter["C"]) >= 1
		),
	];
	return r.some(a => a == true);
}

/**
 * @param {{ A: number; F: number; C: number; N: number; }} counter
 * @param {number} min_caz
 * @returns {boolean}
 */
function strict_CAZyme_definition(counter, min_caz) {
	const r = [
		counter["A"] >= min_caz,
		counter["F"] >= 1,
		counter["C"] >= 1,
	];
	return r.every(a => a == true);
}

/**
 * @param {string} text
 * @param {string} splitter
 * @returns {string[][]}
 */
 function tsv_parse(text, splitter = "\t") {
	let tab = [];
	let tr = [];
	let i = 0;

	text = text.trim().replace(/\r\n/g, "\n");

	while (i < text.length) {
		let c = text[i];
		
		if (c == '"') {
			let val = "";
			
			c = text[++i];//skip "
			while (c != '"' && c != null) {
				val += c;
				c = text[++i];
			}
			c = text[++i];//skip "
			
			tr.push(val.trim());
			
			continue;
		}
		else if (c == splitter) {
			if (text[i - 1] == splitter) {//if this cell is empty
				tr.push("");
			}
			i++;//skip this splitter
			continue;
		}
		else if (c == "\n") {
			tab.push(tr);
			tr = [];
			i++;
			continue;
		}
		else {
			let val = "";
			
			while (c != splitter && c != '\n' && c != null) {
				val += c;
				c = text[++i];
			}
			
			tr.push(val.trim());
			
			continue;
		}
	}
	if (tr.length) {
		tab.push(tr);
	}

	return tab;
}

/**
 * @template T
 * @param {any[][]} table
 * @param {number|string[]} header
 * @param {Partial<{ start_row:number; prototype:object; constructor:any; header_map:(rawHeader:string)=>string; transform: (obj:T,rawHeader:string,val:any)=>T; }>} option
 * @returns {T[]}
 */
function table_to_object_list(table, header, option = { start_row: 0, prototype: null, header_map: null, constructor: null, }) {
	let start_row = option.start_row | 0;
	
	let _table = table.slice(start_row);

	if (header != null) {
		let head;
		if (Number.isSafeInteger(Number(header)) && header < _table.length) {
			for (let i = 0; i < header; ++i) {
				head = _table.shift();
			}
			head = _table.shift();

			// only change name
			if (option.header_map) {
				head = head.map(option.header_map);
			}
		}
		else if (Array.isArray(header)) {
			head = header;
		}
		
		//auto rename duplicate property
		head.reduce((prev, curr, i, arr) => {
			if (prev[arr[i]]) {
				prev[arr[i]]++;
				arr[i] += "_" + prev[arr[i]];
			}
			else {
				prev[arr[i]] = 1;
			}
			return prev;
		}, {});

		//to Object
		if (option.transform) {
			return _table.map((row) => {
				/** @type {T} */
				const obj = option.constructor ? new option.constructor : (option.prototype ? Object.create(option.prototype) : {});
				row.forEach((val, i) => {
					option.transform(obj, head[i], val);
				});
				return obj;
			});
		}
		else {
			return _table.map((row) => {
				return row.reduce((prev, current, i) => {
					prev[head[i]] = current;
					return prev;
				}, option.constructor ? new option.constructor() : (option.prototype ? Object.create(option.prototype) : {}));
			});
		}
	}
	else {
		//@ts-ignore
		return _table;
	}
}

function main() {
	const [
		,// "node",
		,// __filename,
		input,
		output,
	] = process.argv;

	if (fs.existsSync(input) && output != null) {
		CAZyme_cluster(output, input, table_transform, min_caz, max_ns);
	}
	else {
		console.error(`node ${__filename} <input> <ooutput>`)
	}
}

main();

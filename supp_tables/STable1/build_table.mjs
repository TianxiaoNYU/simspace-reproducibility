import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const here = path.dirname(fileURLToPath(import.meta.url));
const workspaceRoot = path.resolve(here, "../../..");
const threadId = "019f9047-2c3a-7852-b086-30ff9814c7d9";
const outputDir = path.join(workspaceRoot, "outputs", threadId);
const renderDir = path.join(outputDir, "STable1_previews");
const outputPath = path.join(outputDir, "Supplementary_Table_1.xlsx");

function parseTsv(text) {
  const normalized = text.replace(/^\uFEFF/, "").replace(/\r\n/g, "\n").trimEnd();
  return normalized.split("\n").map((line) => line.split("\t"));
}

function columnName(index1Based) {
  let value = index1Based;
  let out = "";
  while (value > 0) {
    const remainder = (value - 1) % 26;
    out = String.fromCharCode(65 + remainder) + out;
    value = Math.floor((value - 1) / 26);
  }
  return out;
}

function statusFill(value) {
  if (typeof value !== "string") return null;
  if (value.startsWith("YES")) return "#E7F3EA";
  if (value.startsWith("PARTIAL")) return "#FFF2CC";
  if (value.startsWith("NO")) return "#FCE8E6";
  if (value.startsWith("UNKNOWN")) return "#F4E8FF";
  if (value.startsWith("NOT DESCRIBED")) return "#ECEFF3";
  if (value.startsWith("N/A")) return "#F2F2F2";
  if (value.startsWith("NOT RUNNABLE")) return "#F4CCCC";
  return null;
}

function styleHeader(sheet, rangeAddress) {
  sheet.getRange(rangeAddress).format = {
    fill: "#17324D",
    font: { bold: true, color: "#FFFFFF", size: 10 },
    horizontalAlignment: "center",
    verticalAlignment: "center",
    wrapText: true,
    borders: { preset: "all", style: "thin", color: "#AAB7C4" },
  };
}

function styleBody(sheet, rangeAddress) {
  sheet.getRange(rangeAddress).format = {
    font: { color: "#1F2933", size: 9 },
    verticalAlignment: "top",
    wrapText: true,
    borders: { preset: "all", style: "thin", color: "#D8DEE5" },
  };
}

function addFlatSheet(workbook, name, rows, options = {}) {
  const sheet = workbook.worksheets.add(name);
  sheet.showGridLines = false;
  const rowCount = rows.length;
  const colCount = rows[0].length;
  const endColumn = columnName(colCount);
  const fullAddress = `A1:${endColumn}${rowCount}`;
  sheet.getRange(fullAddress).values = rows;
  styleHeader(sheet, `A1:${endColumn}1`);
  if (rowCount > 1) styleBody(sheet, `A2:${endColumn}${rowCount}`);
  sheet.getRange(`A1:${endColumn}1`).format.rowHeight = options.headerHeight ?? 40;
  if (rowCount > 1) sheet.getRange(`A2:${endColumn}${rowCount}`).format.rowHeight = options.bodyHeight ?? 54;
  sheet.freezePanes.freezeRows(1);
  if (options.freezeColumns) sheet.freezePanes.freezeColumns(options.freezeColumns);
  if (options.tableName) sheet.tables.add(fullAddress, true, options.tableName);
  return { sheet, rowCount, colCount, endColumn };
}

const [mainText, methodsText, evidenceText] = await Promise.all([
  fs.readFile(path.join(here, "STable1.tsv"), "utf8"),
  fs.readFile(path.join(here, "methods.tsv"), "utf8"),
  fs.readFile(path.join(here, "evidence.tsv"), "utf8"),
]);

const mainRows = parseTsv(mainText);
const methodsRows = parseTsv(methodsText);
const evidenceRows = parseTsv(evidenceText);
const expectedMainColumns = 19;

if (!mainRows.every((row) => row.length === expectedMainColumns)) {
  throw new Error(`STable1.tsv must contain exactly ${expectedMainColumns} tab-separated columns on every row.`);
}
if (!methodsRows.every((row) => row.length === methodsRows[0].length)) {
  throw new Error("methods.tsv has inconsistent column counts.");
}
if (!evidenceRows.every((row) => row.length === evidenceRows[0].length)) {
  throw new Error("evidence.tsv has inconsistent column counts.");
}

const evidenceIds = new Set(evidenceRows.slice(1).map((row) => row[0]));
for (const row of mainRows.slice(1)) {
  for (const cell of row.slice(1)) {
    const references = [...cell.matchAll(/E\d{3}/g)].map((match) => match[0]);
    if (references.length === 0) {
      throw new Error(`Every non-label StTable1 cell must cite evidence; missing citation in: ${cell}`);
    }
    for (const evidenceId of references) {
      if (!evidenceIds.has(evidenceId)) {
        throw new Error(`StTable1 cites ${evidenceId}, but that identifier is absent from evidence.tsv.`);
      }
    }
  }
}

await fs.mkdir(renderDir, { recursive: true });

const workbook = Workbook.create();

const main = addFlatSheet(workbook, "STable1", mainRows, {
  freezeColumns: 2,
  headerHeight: 62,
  bodyHeight: 120,
  tableName: "STable1CapabilityAudit",
});

main.sheet.getRange(`A2:B${main.rowCount}`).format = {
  fill: "#EAF2F8",
  font: { bold: true, color: "#17324D", size: 9 },
  verticalAlignment: "top",
  wrapText: true,
  borders: { preset: "all", style: "thin", color: "#D8DEE5" },
};
main.sheet.getRange(`A1:A${main.rowCount}`).format.columnWidth = 18;
main.sheet.getRange(`B1:B${main.rowCount}`).format.columnWidth = 20;
for (const column of ["C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "R", "S"]) {
  main.sheet.getRange(`${column}1:${column}${main.rowCount}`).format.columnWidth = 34;
}
main.sheet.getRange(`Q1:Q${main.rowCount}`).format.columnWidth = 15;

for (let row = 2; row <= main.rowCount; row += 1) {
  for (let col = 3; col <= main.colCount; col += 1) {
    const value = mainRows[row - 1][col - 1];
    const fill = statusFill(value);
    if (fill) main.sheet.getRange(`${columnName(col)}${row}`).format.fill = fill;
  }
}

const methods = addFlatSheet(workbook, "Methods", methodsRows, {
  freezeColumns: 2,
  headerHeight: 42,
  bodyHeight: 66,
  tableName: "STable1Methods",
});
methods.sheet.getRange(`A1:A${methods.rowCount}`).format.columnWidth = 18;
methods.sheet.getRange(`B1:B${methods.rowCount}`).format.columnWidth = 16;
methods.sheet.getRange(`C1:C${methods.rowCount}`).format.columnWidth = 28;
methods.sheet.getRange(`D1:F${methods.rowCount}`).format.columnWidth = 16;
methods.sheet.getRange(`G1:H${methods.rowCount}`).format.columnWidth = 42;
methods.sheet.getRange(`I1:J${methods.rowCount}`).format.columnWidth = 40;
methods.sheet.getRange(`A2:B${methods.rowCount}`).format.fill = "#EAF2F8";
methods.sheet.getRange(`A2:A${methods.rowCount}`).format.font = { bold: true, color: "#17324D" };

const evidence = addFlatSheet(workbook, "Evidence", evidenceRows, {
  freezeColumns: 3,
  headerHeight: 42,
  bodyHeight: 58,
  tableName: "STable1Evidence",
});
const evidenceWidths = [14, 17, 42, 18, 20, 30, 42, 58, 14];
for (let i = 0; i < evidenceWidths.length; i += 1) {
  const col = columnName(i + 1);
  evidence.sheet.getRange(`${col}1:${col}${evidence.rowCount}`).format.columnWidth = evidenceWidths[i];
}
evidence.sheet.getRange(`A2:B${evidence.rowCount}`).format.fill = "#F3F7FA";
evidence.sheet.getRange(`A2:A${evidence.rowCount}`).format.font = { bold: true, color: "#17324D" };

const legend = workbook.worksheets.add("Legend & notes");
legend.showGridLines = false;
legend.getRange("A1:B1").values = [["Supplementary Table 1", "Evidence-audited comparison of spatial simulation capabilities"]];
legend.getRange("A1:B1").format = {
  fill: "#17324D",
  font: { bold: true, color: "#FFFFFF", size: 14 },
  verticalAlignment: "center",
  wrapText: true,
};
legend.getRange("A1:B1").format.rowHeight = 38;
legend.getRange("A3:B10").values = [
  ["Status", "Interpretation"],
  ["YES", "Direct, documented native support in the audited version."],
  ["PARTIAL", "Limited, indirect, optional-backend, or proxy support."],
  ["NO", "The audited source makes the capability incompatible with, or explicitly outside, the method's operating model."],
  ["NOT DESCRIBED", "The capability was not identified in the audited sources; this is not a claim that implementation is impossible."],
  ["UNKNOWN", "Audited sources were incomplete or conflicting."],
  ["N/A", "The capability does not conceptually apply to the method's documented simulation layer."],
  ["NOT RUNNABLE", "Reserved for execution-blocked audits; not used in the current table."],
];
styleHeader(legend, "A3:B3");
styleBody(legend, "A4:B10");
const legendFills = ["#E7F3EA", "#FFF2CC", "#FCE8E6", "#ECEFF3", "#F4E8FF", "#F2F2F2", "#F4CCCC"];
for (let i = 0; i < legendFills.length; i += 1) {
  legend.getRange(`A${i + 4}`).format = {
    fill: legendFills[i],
    font: { bold: true, color: "#17324D" },
    verticalAlignment: "center",
    borders: { preset: "all", style: "thin", color: "#D8DEE5" },
  };
}

legend.getRange("A12:B17").values = [
  ["Audit scope", "Capability/interface audit, not an overall performance ranking."],
  ["Direct spatial molecular effects", "Requires molecular values to vary with coordinate or neighborhood after conceptually holding cell type/base profile fixed; composition-only effects are PARTIAL."],
  ["Cell-cell interaction scope", "Distinguishes planted ligand–receptor/pairwise effects from an explicit receptor-to-downstream-target cascade."],
  ["Truth exports", "Planted latent variables that can be consumed independently by downstream methods."],
  ["Evidence traceability", "Evidence IDs in StTable1 map to versioned primary sources or inspected source code on the Evidence sheet."],
  ["SimSpace release scope", "The SimSpace row audits v0.4.0, which adds a headless CLI and stable CLI outputs without changing the production simulation functions."],
];
styleBody(legend, "A12:B17");
legend.getRange("A12:A17").format = {
  fill: "#DDEBF7",
  font: { bold: true, color: "#17324D" },
  verticalAlignment: "top",
  wrapText: true,
  borders: { preset: "all", style: "thin", color: "#D8DEE5" },
};
legend.getRange("A1:A17").format.columnWidth = 29;
legend.getRange("B1:B17").format.columnWidth = 92;
legend.getRange("A4:B10").format.rowHeight = 36;
legend.getRange("A12:B17").format.rowHeight = 52;
legend.freezePanes.freezeRows(1);

const inspected = {};
const formulaScans = {};
for (const sheetName of ["STable1", "Methods", "Evidence", "Legend & notes"]) {
  inspected[sheetName] = await workbook.inspect({
    kind: "table",
    sheetId: sheetName,
    range: sheetName === "STable1" ? "A1:S9" : undefined,
    include: "values,formulas",
    maxChars: 12000,
  });
  formulaScans[sheetName] = await workbook.inspect({
    kind: "formula",
    sheetId: sheetName,
    maxChars: 2000,
    options: { maxResults: 100 },
  });
}

const formulaErrorPattern = /#REF!|#DIV\/0!|#VALUE!|#NAME\?|#N\/A/;
const formulaScanText = JSON.stringify(formulaScans);
if (formulaErrorPattern.test(formulaScanText)) {
  throw new Error("Formula-error token detected during workbook validation.");
}

for (const sheetName of ["STable1", "Methods", "Evidence", "Legend & notes"]) {
  const preview = await workbook.render({
    sheetName,
    autoCrop: "all",
    scale: sheetName === "STable1" ? 0.65 : 0.8,
    format: "png",
  });
  const safeName = sheetName.replace(/[^A-Za-z0-9]+/g, "_");
  await fs.writeFile(path.join(renderDir, `${safeName}.png`), new Uint8Array(await preview.arrayBuffer()));
}

const xlsx = await SpreadsheetFile.exportXlsx(workbook);
await xlsx.save(outputPath);

const summary = {
  outputPath,
  sheets: ["STable1", "Methods", "Evidence", "Legend & notes"],
  tableRows: mainRows.length - 1,
  evidenceRows: evidenceRows.length - 1,
  mainColumns: mainRows[0].length,
  evidenceCitationsValidated: true,
  formulaErrorScan: "passed; workbook contains no formulas",
  inspected: Object.fromEntries(
    Object.entries(inspected).map(([key, value]) => [key, typeof value === "string" ? value.slice(0, 600) : value]),
  ),
};
console.log(JSON.stringify(summary, null, 2));
